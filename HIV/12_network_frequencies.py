"""Most frequent networks in the squirrel samples, grouped by isomorphism.

Writes 'networks_squirrel/network_frequencies.txt' (pass another directory as
the first argument to write elsewhere). Networks are counted as equal
regardless of how the eNewick is written, using phylozoo's isomorphism tests
(leaf labels are matched; edge lengths and inheritance weights are ignored):

  - rooted samples:   rooted labeled isomorphism  (DirectedPhyNetwork)
  - unrooted samples: semi-directed isomorphism   (SemiDirectedPhyNetwork)

Each network is compared to one representative per class found so far and
either joins that class or starts a new one. Classes are ranked by size (ties:
first seen) and the TOP largest are listed with their first member's eNewick.

Requires (version used): Python 3.11, phylozoo 0.2.6 (pip install phylozoo==0.2.6).
"""

import sys
from pathlib import Path

from phylozoo import DirectedPhyNetwork, SemiDirectedPhyNetwork
from phylozoo.core.network.dnetwork.isomorphism import is_isomorphic as d_iso
from phylozoo.core.network.sdnetwork.isomorphism import is_isomorphic as sd_iso

DIR = Path(__file__).parent / "networks_squirrel"
OUTFILE = Path(sys.argv[1] if len(sys.argv) > 1 else DIR) / "network_frequencies.txt"
TOP = 5
# (title, sample file, network class, isomorphism test)
RUNS = [
    ("Rooted, all subtypes",      "squirrel_rooted_sampling_HIV_all_subtypename.txt",     DirectedPhyNetwork,     d_iso),
    ("Rooted, without J and H",   "squirrel_rooted_sampling_HIV_noJH_subtypename.txt",    DirectedPhyNetwork,     d_iso),
    ("Unrooted, all subtypes",    "squirrel_unrooted_sampling_HIV_all_subtypename.txt",   SemiDirectedPhyNetwork, sd_iso),
    ("Unrooted, without J and H", "squirrel_unrooted_sampling_HIV_noJH_subtypename.txt",  SemiDirectedPhyNetwork, sd_iso),
]

HEADER = """\
Five most frequent networks per squirrel sampling run, for the rooted and the
unrooted sampling.

Rooted runs are grouped by rooted labeled isomorphism (phylozoo 0.2.6
DirectedPhyNetwork). Unrooted runs are grouped by semi-directed isomorphism
(phylozoo 0.2.6 SemiDirectedPhyNetwork). In both cases networks are counted as
equal regardless of how the eNewick is written.
"""


def group(lines, cls, iso):
    """Return [[count, first eNewick, representative network], ...], largest first."""
    classes = []
    for line in lines:
        net = cls.from_string(line, format="enewick")
        for c in classes:
            if iso(net, c[2]):
                c[0] += 1
                break
        else:
            classes.append([1, line, net])
    return sorted(classes, key=lambda c: -c[0])     # stable: ties keep first-seen order


out = [HEADER]
for title, file, cls, iso in RUNS:
    lines = [l.strip() for l in (DIR / file).read_text().splitlines() if l.strip()]
    n = len(lines)
    classes = group(lines, cls, iso)
    top = sum(c[0] for c in classes[:TOP])
    out.append(f"\n{'=' * 72}\n{title}\n{'=' * 72}")
    out.append(f"Input:    {file}\nNetworks: {n}\nDistinct: {len(classes)}")
    out.append(f"Top {TOP}:    {top}/{n} ({100 * top / n:.1f}%)\n")
    for i, (count, enewick, _) in enumerate(classes[:TOP], 1):
        out.append(f"#{i}  {count}/{n}  ({100 * count / n:.1f}%)\n{enewick}\n")

OUTFILE.write_text("\n".join(out))
print("wrote", OUTFILE)
