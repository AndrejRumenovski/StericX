"""Render exact self instruction counts from Callgrind --separate-callers=40.

Widths are modeled instruction counts, not CPU time. Callgrind encodes its
caller context as apostrophe-separated function names. Inline instructions
remain charged to their containing function. No edge costs are added twice.
"""

import argparse
import collections
import hashlib
import html
import json
import re
from pathlib import Path


def parse(path):
    names, costs = {}, collections.Counter()
    current, call, positions, total = None, False, 1, None
    for line in path.read_text().splitlines():
        if line.startswith("positions:"):
            positions = len(line.split()) - 1
        elif line.startswith("events:"):
            assert line.split()[1:] == ["Ir"], "one instruction event required"
        elif line.startswith(("fn=", "cfn=")):
            prefix, value = line.split("=", 1)
            match = re.fullmatch(r"\((\d+)\)(?: (.*))?", value)
            if match:
                identity, name = match.groups()
                if name is not None:
                    names[identity] = name
                value = identity
            else:
                names[value] = value
            if prefix == "fn":
                current, call = value, False
        elif line.startswith("calls="):
            call = True
        elif line.startswith("totals:"):
            total = int(line.split()[1])
        elif line and line[0] in "0123456789+-*" and current is not None:
            fields = line.split()
            if not call:
                costs[current] += int(fields[positions])
            call = False
    assert total is not None and sum(costs.values()) == total, "self costs must reconstruct total"
    return {names[k]: value for k, value in costs.items() if value}, total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    costs, total = parse(args.input)
    root = {"name": "all instructions", "value": total, "children": {}}
    for context, value in costs.items():
        node = root
        for name in reversed(context.split("'")):
            node = node["children"].setdefault(name, {"name": name, "value": 0, "children": {}})
            node["value"] += value

    def depth(node):
        return 1 + max((depth(c) for c in node["children"].values()), default=0)

    height = 100 + 20 * depth(root)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="{height}" viewBox="0 0 1400 {height}">',
        '<rect width="100%" height="100%" fill="#faf9f6"/>',
        '<text x="20" y="25" font-family="sans-serif" font-size="18">StericX: Callgrind instruction flamegraph</text>',
        f'<text x="20" y="48" font-family="sans-serif" font-size="12">{total:,} modeled instructions; caller contexts limited to 40 frames; widths are not CPU time. Hover for counts.</text>']

    def draw(node, x, level):
        width = 1360 * node["value"] / total
        y = height - 25 - 20 * level
        color = int(hashlib.sha256(node["name"].encode()).hexdigest()[:2], 16)
        if width >= 0.15:
            name = html.escape(node["name"])
            parts.append(f'<g><title>{name}: {node["value"]:,} instructions ({node["value"]/total:.3%})</title>'
                f'<rect x="{x:.3f}" y="{y}" width="{width:.3f}" height="18" fill="rgb(245,{110+color//3},{60+color//5})" stroke="#fff" stroke-width="0.3"/>')
            if width > 36:
                label = html.escape(node["name"][:max(1, int(width / 6.4) - 2)])
                parts.append(f'<text x="{x+3:.3f}" y="{y+13}" font-family="monospace" font-size="10">{label}</text>')
            parts.append('</g>')
        for child in sorted(node["children"].values(), key=lambda n: n["name"]):
            draw(child, x, level+1)
            x += 1360 * child["value"] / total

    draw(root, 20, 0)
    parts.append('</svg>')
    args.output.write_text("\n".join(parts) + "\n")
    args.output.with_suffix('.json').write_text(json.dumps({
        "input": str(args.input), "sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "total_instructions": total, "self_costs": costs,
        "semantics": __doc__, "depth": depth(root)}, indent=2) + "\n")
    print(json.dumps({"total_instructions": total, "contexts": len(costs), "depth": depth(root)}))


if __name__ == '__main__':
    main()
