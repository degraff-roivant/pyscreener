from pathlib import Path

import pytest

from pyscreener.docking.smina import SminaRunner

TEST_DIR = Path(__file__).parents[0]
TEST_DATA_DIR = TEST_DIR / "data"


@pytest.fixture
def multi_log():
    return TEST_DATA_DIR / "smina_multi_log.txt"


@pytest.mark.parametrize(
    "contents,scores_true",
    [
        (
            r"""Weights      Terms
-0.035579    gauss(o=0,_w=0.5,_c=8)
-0.005156    gauss(o=3,_w=2,_c=8)
0.840245     repulsion(o=0,_c=8)
-0.035069    hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439    non_dir_h_bond(g=-0.7,_b=0,_c=8)
1.923        num_tors_div

Using random seed: 479703264""",
            None,
        ),
        (
            r"""-0.035069    hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439    non_dir_h_bond(g=-0.7,_b=0,_c=8)
1.923        num_tors_div

Using random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -4.3       0.000      0.000    
2       -4.3       0.084      2.626    """,
            [-4.3, -4.3],
        ),
        (
            r"""ing random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -4.3       0.000      0.000    
2       -4.8       0.084      2.626    
7       -4.0       1.726      2.810    
8       -3.9       9.133      10.292   
9       -3.4       9.149      10.463   
Using random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -3.2       0.000      0.000    
2       -3.2       0.029      1.837    
3       -3.1       1.226      1.487 """,
            [-4.3, -4.8, -4.0, -3.9, -3.4],
        ),
    ],
)
def test_parse_log(tmp_path: Path, contents, scores_true):
    logfile = tmp_path / "log.txt"
    logfile.write_text(contents)

    scores_actual = SminaRunner.parse_logfile(logfile)

    assert scores_actual == pytest.approx(scores_true)


@pytest.mark.parametrize(
    "contents,scoress_true",
    [
        (
            r"""   _______  _______ _________ _        _______ 
  (  ____ \(       )\__   __/( (    /|(  ___  )
  | (    \/| () () |   ) (   |  \  ( || (   ) |
  | (_____ | || || |   | |   |   \ | || (___) |
  (_____  )| |(_)| |   | |   | (\ \) ||  ___  |
        ) || |   | |   | |   | | \   || (   ) |
  /\____) || )   ( |___) (___| )  \  || )   ( |
  \_______)|/     \|\_______/|/    )_)|/     \|


smina is based off AutoDock Vina. Please cite appropriately.

Weights      Terms
-0.035579    gauss(o=0,_w=0.5,_c=8)
-0.005156    gauss(o=3,_w=2,_c=8)
0.840245     repulsion(o=0,_c=8)
-0.035069    hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439    non_dir_h_bond(g=-0.7,_b=0,_c=8)
1.923        num_tors_div

Using random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -4.3       0.000      0.000    
2       -4.3       0.084      2.626    
5       -4.3       0.119      2.123    
6       -4.3       0.134      2.122    
7       -4.2       1.726      2.810    
9       -4.2       9.149      10.463   
Using random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -3.2       0.000      0.000    
3       -3.1       1.226      1.487    
4       -3.1       11.736     11.911   
5       -3.0       5.558      5.863    
6       -3.0       3.729      4.105    
7       -2.9       12.243     12.418   
8       -2.8       3.963      4.093    
9       -2.8       12.655     12.769   
Using random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -3.2       0.000      0.000    
4       -3.1       8.733      9.740    
5       -3.0       8.631      9.962    
6       -3.0       10.735     11.360   
8       -2.8       10.824     11.502   
9       -2.7       1.312      2.526    
Using random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -7.6       0.000      0.000    
2       -7.3       2.023      3.089    
6       -6.7       2.955      4.650    
7       -6.6       1.994      2.468    
9       -6.5       4.548      8.903    
Using random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -2.8       0.000      0.000    
2       -2.7       1.242      1.242    
7       -2.5       10.708     11.250   
8       -2.5       10.376     10.663   
9       -2.4       2.528      3.096    
Using random seed: 479703264

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -5.9       0.000      0.000    
2       -5.9       6.311      8.222    
7       -5.6       7.803      9.431    
8       -5.4       1.394      2.181    
9       -5.4       1.658      2.659    
""",
            [
                [-4.3, -4.3, -4.3, -4.3, -4.2, -4.2],
                [-3.2, -3.1, -3.1, -3.0, -3.0, -2.9, -2.8, -2.8],
                [-3.2, -3.1, -3.0, -3.0, -2.8, -2.7],
                [-7.6, -7.3, -6.7, -6.6, -6.5],
                [-2.8, -2.7, -2.5, -2.5, -2.4],
                [-5.9, -5.9, -5.6, -5.4, -5.4],
            ],
        )
    ],
)
def test_batch_parse_log(tmp_path: Path, contents, scoress_true):
    logfile = tmp_path / "log.txt"
    logfile.write_text(contents)

    scoress_actual = SminaRunner.batch_parse_logfile(logfile)

    for scores_actual, scores_true in zip(scoress_actual, scoress_true):
        assert scores_actual == pytest.approx(scores_true)
