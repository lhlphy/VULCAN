#!/usr/bin/env python3
"""Run one KIDA species search and save source-specific raw HTML fragments."""

from __future__ import annotations

import argparse
import html
import json
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://kida.astrochem-tools.org/"
USER_AGENT = "VULCAN-Reaction-Search/1.0 manual-research"
ENDPOINTS = {
    "unimolecular": "searchCrReactions",
    "bimolecular": "searchBimoReactions",
    "termolecular": "searchTermoReactions",
    "surface": "searchSurfaceReactions",
}


def default_temp_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "data_temp"


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.+-]+", "_", value).strip("_") or "reaction"


def request_bytes(url: str, data: bytes | None = None) -> bytes:
    request = urllib.request.Request(
        url,
        data=data,
        headers={"User-Agent": USER_AGENT, "X-Requested-With": "XMLHttpRequest"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def extract_species_map(page: str) -> dict[str, object]:
    match = re.search(r"'tabSpeciesIsotopes':\s*'(\{.*?\})'", page)
    if not match:
        raise RuntimeError("KIDA species mapping was not found; inspect the saved search HTML")
    return json.loads(html.unescape(match.group(1)))


def text_preview(fragment: str) -> str:
    text = re.sub(r"<[^>]+>", " ", fragment)
    text = html.unescape(text).replace("\u2192", "->")
    return " ".join(text.split())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, help='KIDA query such as "OH + CH4"')
    parser.add_argument("--reactprod", choices=("reactants", "products", "both"), default="reactants")
    parser.add_argument("--domain", choices=("Astro", "Planeto", "Both"), default="Both")
    parser.add_argument("--ionneutral", choices=("ion", "neutral"), default="neutral")
    parser.add_argument("--output-dir", type=Path, default=default_temp_dir())
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    query = urllib.parse.urlencode(
        {
            "kida[species]": args.species,
            "kida[reactprod]": args.reactprod,
            "kida[ionneutral]": args.ionneutral,
            "kida[astroplaneto]": args.domain,
        }
    )
    search_url = urllib.parse.urljoin(BASE_URL, "search") + "?" + query
    search_page = request_bytes(search_url).decode("utf-8", errors="replace")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    stem = f"kida_{slug(args.species)}_{stamp}"
    initial_path = args.output_dir / f"{stem}_search.html"
    initial_path.write_text(search_page, encoding="utf-8")

    species_map = extract_species_map(search_page)
    post = urllib.parse.urlencode(
        {
            "tabSpeciesIsotopes": json.dumps(species_map, separators=(",", ":")),
            "reactprod": args.reactprod,
            "ionneutral": args.ionneutral,
            "astroplaneto": args.domain,
            "query": args.species,
            "isomers": "1",
            "ids": "",
            "page": "1",
        }
    ).encode()

    print(f"Saved KIDA search page: {initial_path}")
    for label, endpoint in ENDPOINTS.items():
        fragment = request_bytes(urllib.parse.urljoin(BASE_URL, endpoint), post).decode(
            "utf-8", errors="replace"
        )
        path = args.output_dir / f"{stem}_{label}.html"
        path.write_text(fragment, encoding="utf-8")
        preview = text_preview(fragment)
        status = preview[:300] if preview else "no results"
        print(f"{label}: {status}")
        print(f"  raw: {path}")

    print("Inspect KIDA detail pages and expert/method labels before selecting a value.")


if __name__ == "__main__":
    main()
