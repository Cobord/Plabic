"""
test for sign patterns
"""

import random
from typing import List
from plabic.sign_pattern import SignPattern

def test_degenerate_sign_pattern():
    """
    too few nonzeros
    """
    errored = False
    t = []
    try:
        t_made = SignPattern(t)
    except ValueError:
        errored = True
    assert errored

    errored = False
    t = [0]
    try:
        t_made = SignPattern(t)
    except ValueError:
        errored = True
    assert errored

    errored = False
    t = [0,0]
    try:
        t_made = SignPattern(t)
    except ValueError:
        errored = True
    assert errored

    errored = False
    t = [0,1]
    try:
        t_made = SignPattern(t)
    except ValueError:
        errored = True
    assert not errored
    assert t_made.var == 0
    assert t_made.varbar == 1

def test_typical():
    """
    something that isn't an edge case
    """
    t = [4,-1,0,-2]
    t_made = SignPattern(t)
    assert t_made.var == 1
    assert t_made.varbar == 3

def test_alt():
    """
    property based test
    that var(v) + varbar(alt(v)) = n-1
    for all v in R^n that is not all 0
    """
    my_n = 4
    trials_num = 30
    word_lengths = [0,1,2,3,5,10]
    for _ in range(trials_num):
        random_word : List[float] = \
            random.choices(range(-my_n,my_n+1), k = random.choice(word_lengths))
        if random_word == [0]*len(random_word):
            continue
        random_sign_pattern = SignPattern(random_word)
        alt_random_word = [((-1)**i)*v for (i,v) in enumerate(random_word)]
        random_sign_pattern_alt = SignPattern(alt_random_word)
        assert random_sign_pattern.var + random_sign_pattern_alt.varbar == \
            len(random_word)-1, f"{random_word} {alt_random_word}"
