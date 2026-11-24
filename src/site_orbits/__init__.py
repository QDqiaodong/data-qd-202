"""Domain data processing API; synthetic offline baseline."""
from .service import process
from .supercell import assemble_supercell_shards, build_supercell_shard
from .extract import extract_supercell_groups
from .corpus import Corpus, CorpusError, FoldConflict
from .draft import solve_draft
from .handoff import handoff_windows, rebuild_handoff
