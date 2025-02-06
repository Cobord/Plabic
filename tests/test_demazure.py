"""
test for hopping operators and Demazure products
"""

from plabic import OneLinePermutation, hopping_operator, up_arrow

def test_example():
    """
    test the example provided in paper
    """
    w : OneLinePermutation = [8,9,1,7,2,6,4,3,5]
    w_final = hopping_operator(1, list(range(8,1,-1)), w)
    assert w_final == [8,9,7,6,2,5,4,3,1]

    w : OneLinePermutation = [8,9,1,7,2,6,4,3,5]
    # [8,9,1,7,2,6,4,3,5]
    # [8,9,7,1,2,6,4,3,5]
    # [8,9,7,6,2,1,4,3,5]
    # [8,9,7,6,2,5,4,3,1]
    w_final = hopping_operator(1, [2,7,5,6,3], w)
    assert w_final == [8,9,2,7,5,6,4,3,1]

def test_up_operator():
    """
    test the example provided in paper
    """
    w : OneLinePermutation = [8,9,1,7,2,6,4,3,5]
    l_answer = up_arrow(w, 2)
    assert l_answer == [7,9,8]

    w : OneLinePermutation = [8,9,1,7,2,6,4,3,5]
    l_answer = up_arrow(w, 4)
    assert l_answer == [6,7,9,8]
