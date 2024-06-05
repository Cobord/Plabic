"""
test for cyclic utilities
"""
#pylint:disable=unused-import
from plabic import generate_lyndons

def test_unary_lyndon():
    """
    lyndon words in one letters
    """
    alphabet = ["0"]
    expected_lyndons = ["0"]
    all_observed_lyndons = ["".join(word) for word in generate_lyndons(alphabet,5)]
    assert all_observed_lyndons == expected_lyndons

def test_binary_lyndon():
    """
    lyndon words in two letters
    """
    alphabet = ["0","1"]
    expected_lyndons = ["0","1","01","001","011","0001","0011","0111",
                        "00001", "00011", "00101", "00111", "01011", "01111"]
    for idx,word in enumerate(generate_lyndons(alphabet,5)):
        real_word = "".join(word)
        assert real_word == expected_lyndons[idx]
        if idx == len(expected_lyndons):
            break
