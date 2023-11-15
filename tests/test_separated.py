"""
test strongly/weakly separated collections of sets
"""
from plabic.separated_sets import StronglySeparatedCollection,WeaklySeparatedCollection

def test_strong_disproof():
    """
    test of disproof of proposition that two sets are strongly separated
    """
    disproof = StronglySeparatedCollection.strongly_separated_disproof(set([1,3]),set([2]))
    assert disproof is not None
    assert disproof == (1,2,3)

def test_strongly_separated_sets():
    """
    test of strongly separated collections of sets
    """
    three_ex = StronglySeparatedCollection(
        [set(),set([1]),set([1,2]),
         set([1,2,3]),
         set([3]),set([3,2])],my_min=1,my_max=3)
    assert not three_ex.maximal_by_size()
    assert three_ex.my_min == 1
    assert three_ex.my_max == 3
    three_ex.append(set([2]))
    # pylint:disable=invalid-name
    num_sets = len(three_ex)
    three_ex.append(set([2]))
    three_ex.append(set([2]))
    assert len(three_ex) == num_sets, f"{num_sets} vs {len(three_ex)}"
    assert three_ex.maximal_by_size()
    try:
        three_ex.append(set([5,4]))
        errored = False
    except ValueError as e:
        errored = True
        #pylint:disable=line-too-long
        exp_msg = "{4, 5} was not a subset of [1, 3]"
        assert str(e) == exp_msg, f"{str(e)} vs {exp_msg}"
    assert errored

def test_weak_disproof():
    """
    test of disproof of proposition that two sets are weakly separated
    """
    disproof_two = WeaklySeparatedCollection.weakly_separated_disproof(set([1,3]),set([2]))
    assert disproof_two is None
    disproof_two = WeaklySeparatedCollection.weakly_separated_disproof(set([1,3]),set([2,4]))
    assert disproof_two is not None
    assert disproof_two == (1,2,3,4)

def test_nonuniform_weakly_separated_sets():
    """
    test of non-uniform weakly separated collections of sets
    """
    four_ex = WeaklySeparatedCollection(
        [set(),set([1]),set([1,2]),
         set([1,2,3]),
         set([1,2,3,4]),
         set([4]),set([4,3]),
         set([4,3,2]),
         set([2,3]),set([2]),set([3]),
         set([4,1]),set([3,4,1]),set([4,1,2])],uniform_k=False)
    assert not four_ex.maximal_by_size()
    four_ex.append(set([1,3]))
    assert four_ex.my_min == 1
    assert four_ex.my_max == 4
    assert four_ex.all_sets_k is None, f"{four_ex.all_sets_k}"
    assert four_ex.maximal_by_size()
    # pylint:disable=invalid-name
    num_sets = len(four_ex)
    four_ex.append(set([2]))
    four_ex.append(set([2]))
    assert len(four_ex) == num_sets, f"{num_sets} vs {len(four_ex)}"

def test_uniform_weakly_separated_sets():
    """
    test of uniform weakly separated collections of sets
    """
    four_ex = WeaklySeparatedCollection(
        [set([1,2]),
         set([4,3]),
         set([2,3]),
         set([4,1])],uniform_k=True)
    assert four_ex.my_min == 1
    assert four_ex.my_max == 4
    assert four_ex.all_sets_k == 2, f"{four_ex.all_sets_k}"
    assert not four_ex.maximal_by_size()
    four_ex.append(set([1,3]))
    assert four_ex.maximal_by_size()
    # pylint:disable=invalid-name
    num_sets = len(four_ex)
    try:
        four_ex.append(set([2,4]))
        errored = False
    except ValueError as e:
        errored = True
        #pylint:disable=line-too-long
        exp_msg = "{2, 4} was not weakly separated with something because we are already at maximum 5"
        assert str(e) == exp_msg, f"{str(e)} vs {exp_msg}"
    assert errored
    four_ex.append(set([1,3]))
    assert len(four_ex) == num_sets, f"{num_sets} vs {len(four_ex)}"
