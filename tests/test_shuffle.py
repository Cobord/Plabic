"""
test for shuffles
"""

from plabic.shuffle import inverse_of_shuffles, inverse_perm, is_inverse_shuffle, shuffle_product

def test_is_shuffle():
    """
    test if something is or not a inverse of a shuffle
    """
    trying = [0,1,2,3,4]
    assert is_inverse_shuffle(trying,2,3)
    trying = [0,2,1,3,4]
    assert is_inverse_shuffle(trying,2,3)

def test_shuffle_gen():
    """
    test generation of shuffles
    """
    expected_shuffles = ["01234","02134","02314","02341",\
                         "23401","23041","23014","20341","20134","20314"]
    shuffle_count = 0
    for sigma_inverse in inverse_of_shuffles(2,3):
        sigma = inverse_perm(sigma_inverse, 5)
        sigma_str = "".join((str(x) for x in sigma))
        assert sigma_str in expected_shuffles
        shuffle_count += 1
    assert shuffle_count == len(expected_shuffles)

def test_example():
    """
    basic example of shuffle product
    """
    factor_1 = ["a","b"]
    factor_2 = ["x","y"]
    expected_shuffles = ["abxy","axby","xaby","axyb","xayb","xyab"]
    summand_count = 0
    for summand in shuffle_product(factor_1, factor_2):
        summand_str = "".join(summand)
        assert summand_str in expected_shuffles
        summand_count += 1
    assert summand_count == len(expected_shuffles)

def test_repeated_summand():
    """
    same summand repeated
    """
    factor_1 = ["a","a","a"]
    factor_2 = ["a","a"]
    expected_shuffles = ["a"*(len(factor_1) + len(factor_2))]
    summand_count = 0
    for summand in shuffle_product(factor_1, factor_2):
        summand_str = "".join(summand)
        assert summand_str in expected_shuffles
        summand_count += 1
    assert summand_count == 10

def test_self_shuffle():
    """
    one object shuffled with itself so
    it needs cloning first
    """
    factor = ["a","a","a"]
    expected_shuffles = ["a"*(len(factor) + len(factor))]
    summand_count = 0
    for summand in shuffle_product(factor, factor):
        summand_str = "".join(summand)
        assert summand_str in expected_shuffles
        summand_count += 1
    assert summand_count == 20
