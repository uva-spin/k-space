"""Explicit resource options for the historical physical demonstration clients."""
from pathlib import Path

def add_pdf_arguments(parser):
    parser.add_argument('--library', type=Path, help='Path to the compiled libkt_lhapdf bridge')
    parser.add_argument('--pdf-data', type=Path, help='Directory containing the pinned LHAPDF set')

def pdf_paths(args, repo):
    # Retain the historical defaults for existing local reproduction scripts.
    # Public recipes always pass explicit paths and have no checkout-name dependency.
    runs = repo.parent.parent
    library = args.library or runs/'m2_prescription_2026-09-30/runtime/libkt_lhapdf.dylib'
    data = args.pdf_data or runs/'m1_resolution_2026-09-30/runtime/lhapdf/share/LHAPDF'
    if not library.is_file() or not data.is_dir():
        raise RuntimeError('PDF runtime not found; supply --library and --pdf-data explicitly')
    return library.resolve(), data.resolve()
