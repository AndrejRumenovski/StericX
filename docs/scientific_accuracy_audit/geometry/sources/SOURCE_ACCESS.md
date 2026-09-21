# Geometry source access and limits

The original HTTP captures from 19 September 2026 remain unchanged in
[manifest.json](manifest.json), including unsuccessful responses. An HTTP status
of 200 alone does not establish access to an original article's full text.

| Reference | Evidence actually available | Limit |
|---|---|---|
| Verloop et al. (1976), DOI 10.1016/B978-0-12-060307-7.50010-9 | Bibliographic record; official Morfeus/AaronTools conventions; frozen runnable Morfeus source | Original chapter text and historical radius parameterization were not independently inspected. |
| Radhakrishnan and Agranat (1991), DOI 10.1007/BF00676621 | Publisher page; official Morfeus equations/source; independent analytic geometry | Original full-text derivation was unavailable. Agreement establishes the runnable convention, not definitive historical textual reproduction. |
| Bondi (1964), DOI 10.1021/j100785a001 | Bibliographic attribution; independently versioned software radius tables | Publisher and archived-paper retrieval failed. Br 1.85 versus 1.83 Å remains a table mismatch without direct adjudication against the original printed table. |
| Cordero et al. (2008), DOI 10.1039/B801115J | Publisher abstract available through indexed [RSC article record](https://pubs.rsc.org/en/content/articlepdf/2008/dt/b801115j), checked 21 September 2026 | The frozen HTML request received 403 and the attempted supporting-data URL 404. The abstract supports crystallographic derivation of radii, not a universal 1.3-times binary bonding cutoff. |
| Mantina et al. (2009), DOI 10.1021/jp8111556 | Indexed [ACS abstract](https://pubs.acs.org/doi/abs/10.1021/jp8111556), checked 21 September 2026, explicitly gives the extension B = 1.92 Å | Original PMC request failed; no full-text table reconstruction is claimed. This primary abstract supports the boron extension, not every StericX radius. |
| SambVca official manual | Full frozen manual HTTP body | Descriptor defaults are reference conventions, not a uniquely correct atomic-surface model. |
| Morfeus 0.8.0 | Frozen installed source, official documentation and separately executed calculations | Independent implementation, not physical truth. Its B5 angular approximation and odd-grid nonadditivity are retained. |
| glam 0.30.10 | Exact installed SSE2 quaternion source, hash and API documentation in [alignment/](../alignment/) | This is a dependency of the system under investigation, not an independent scientific implementation. |

Only short factual paraphrases of publisher abstracts are used here. Source
access limits remain explicit in the claim confidence and do not become passes
because a software test agrees.
