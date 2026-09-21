"""Generate the final scope-correct model claim inventory.

The first draft is retained in models/results/claims_draft_preserved.json;
models_polish_claims applies the documented evidence corrections and wording.
"""
from models_polish_claims import main

if __name__ == "__main__":
    main()
