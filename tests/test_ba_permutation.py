"""
test for (bounded) affine permutations
"""
import itertools
from math import factorial
from plabic import AffinePermutation, BoundedAffinePermutation

def test_coxeter_muls():
    """
    test coxeter
    """
    trial_n = 5
    cox_gens = {}
    for idx in range(trial_n):
        s_idx = AffinePermutation(n_val=trial_n,coxeter_word=[idx])
        if idx==0:
            assert all(s_idx[jdx]==jdx for jdx in range(2,trial_n))
            assert s_idx[1]==0
            assert s_idx[trial_n]==trial_n+1
        else:
            assert all(jdx in [s_idx[jdx],idx,idx+1] for jdx in range(1,trial_n+1))
            assert s_idx[idx]==idx+1
            assert s_idx[idx+1]==idx
        cox_gens[idx] = s_idx
    for (idx,jdx) in itertools.permutations(range(trial_n),2):
        s_idx = cox_gens[idx]
        s_jdx = cox_gens[jdx]
        combined_obs = AffinePermutation(n_val=trial_n,coxeter_word=[idx,jdx])
        combined_exp = s_idx*s_jdx
        assert combined_exp.my_coxeter_word == [idx,jdx]
        assert combined_obs.my_coxeter_word == [idx,jdx]
        assert combined_exp == combined_obs
    for (idx,jdx,kdx) in itertools.permutations(range(trial_n),3):
        s_idx = cox_gens[idx]
        s_jdx = cox_gens[jdx]
        s_kdx = cox_gens[kdx]
        combined_obs = AffinePermutation(n_val=trial_n,coxeter_word=[idx,jdx,kdx])
        combined_exp = s_idx*s_jdx*s_kdx
        assert combined_exp.my_coxeter_word == [idx,jdx,kdx]
        assert combined_obs.my_coxeter_word == [idx,jdx,kdx]
        assert combined_exp == combined_obs

def test_bounded_affine():
    """
    construct a bounded affine permutation
    """
    _my_ba = BoundedAffinePermutation(some_vals = {1:3,2:2})

def test_eq_coxeter():
    """
    equality comparison giving true
    makes sure both have the shortest coxeter word that
    is a presentation for that element
    """
    trial_n = 5
    cox_gens = {}
    for idx in range(trial_n):
        s_idx = AffinePermutation(n_val=trial_n,coxeter_word=[idx])
        if idx==0:
            assert all(s_idx[jdx]==jdx for jdx in range(2,trial_n))
            assert s_idx[1]==0
            assert s_idx[trial_n]==trial_n+1
        else:
            assert all(jdx in [s_idx[jdx],idx,idx+1] for jdx in range(1,trial_n+1))
            assert s_idx[idx]==idx+1
            assert s_idx[idx+1]==idx
        cox_gens[idx] = s_idx
    actually_s0 = AffinePermutation(some_vals = {1:0,2:2,3:3,4:4,5:6})
    assert actually_s0.my_coxeter_word is None
    assert actually_s0 == cox_gens[0]
    assert actually_s0.my_coxeter_word == [0]

def test_all_bounded_affine():
    """
    Get all sum_k=0^n n!/k! bounded affine permutations
    when do the full iterator
    """
    expected_counts = [2,5,16,65,326,1957]
    for cur_n,expected_count in enumerate(expected_counts,1):
        count = 0
        for _my_ba in BoundedAffinePermutation.all_bounded_affine_perms(cur_n):
            count += 1
        assert count == expected_count

def test_all_affine_perms():
    """
    Get all affine permutations up to given length
    when do the full iterator
    """
    # coefficient of q^1 in (1-q^n)/(1-q)^n as n goes from 2 to 7
    expected_counts_1 = [2,3,4,5,6,7]
    # coefficient of q^2 in (1-q^n)/(1-q)^n as n goes from 2 to 7
    expected_counts_2 = [2,6,10,15,21,28]
    for cur_n,(expected_count_1,expected_count_2) in \
        enumerate(zip(expected_counts_1,expected_counts_2),start=2):
        count = 0
        for len_my_ba,my_ba in AffinePermutation.all_affine_perms(cur_n,0):
            count += 1
            assert my_ba.is_lift_from_sn
            assert len(my_ba)<=0
            assert len_my_ba==len(my_ba)
        assert count == 1
        count = 0
        for len_my_ba,my_ba in AffinePermutation.all_affine_perms(cur_n,1):
            count += 1
            assert len(my_ba)<=1
            assert len_my_ba==len(my_ba)
        assert count == expected_count_1+1
        count = 0
        for len_my_ba,my_ba in AffinePermutation.all_affine_perms(cur_n,2):
            count += 1
            assert len(my_ba)<=2
            assert len_my_ba==len(my_ba)
        assert count == expected_count_2+expected_count_1+1, f"{cur_n} {count}"

def test_all_finite_perms():
    """
    Get all finite permutations up to given length
    when do the full iterator
    """
    # number of permutations with 1 inversion
    expected_counts_1 = [1,2,3,4,5]
    # number of permutations with 2 inversions
    expected_counts_2 = [0,2,5,9,14]
    for cur_n,(expected_count_1,expected_count_2) in \
        enumerate(zip(expected_counts_1,expected_counts_2),start=2):
        count = 0
        for len_my_ba,my_ba in AffinePermutation.all_finite_perms(cur_n,0):
            count += 1
            assert my_ba.is_lift_from_sn
            assert len(my_ba)<=0
            assert len_my_ba==len(my_ba)
        assert count == 1
        count = 0
        for len_my_ba,my_ba in AffinePermutation.all_finite_perms(cur_n,1):
            count += 1
            assert my_ba.is_lift_from_sn
            assert len(my_ba)<=1
            assert len_my_ba==len(my_ba)
        assert count == expected_count_1+1
        count = 0
        for len_my_ba,my_ba in AffinePermutation.all_finite_perms(cur_n,2):
            count += 1
            assert my_ba.is_lift_from_sn
            assert len(my_ba)<=2
            assert len_my_ba==len(my_ba)
        assert count == expected_count_2+expected_count_1+1
        count = 0
        longest_count = [0,0]
        for len_my_ba,my_ba in AffinePermutation.all_finite_perms(cur_n):
            count += 1
            assert my_ba.is_lift_from_sn
            assert len_my_ba==len(my_ba)
            if len_my_ba>longest_count[0]:
                longest_count[0] = len_my_ba
                longest_count[1] = 1
            elif len_my_ba==longest_count[0]:
                longest_count[1] += 1
        assert count == factorial(cur_n)
        assert longest_count[0] == cur_n*(cur_n-1)//2
        assert longest_count[1] == 1

def test_all_parabolic_perms():
    """
    Get all parabolic permutations
    when do the full iterator
    """
    my_n = 9
    my_k = 5
    # all 120 of S_{my_k = 5} by length
    sk_by_length = [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1]
    # all 24 of S_{my_n-my_k = 4} by length
    s_nmk_by_length = [1, 3, 5, 6, 5, 3, 1]
    max_len = len(sk_by_length)+len(s_nmk_by_length)-2
    counts = {z_var:0 for z_var in range(max_len+1)}
    for len_my_ba,my_ba in AffinePermutation.all_k_parabolic_perms(my_n,my_k,max_len):
        counts[len_my_ba] += 1
        assert my_ba.is_lift_from_sn
        assert len(my_ba)<=max_len
        assert len_my_ba==len(my_ba)
        is_grassmannian, which_k = my_ba.is_k_grassmannian()
        assert not is_grassmannian or which_k!=my_k
        all_descents = my_ba.right_descents()
        assert my_k not in all_descents
    for my_len in range(max_len+1):
        kdx_max = min(my_len,len(sk_by_length)-1)
        kdx_min = max(my_len-len(s_nmk_by_length)+1,0)
        expected_count = \
            sum(sk_by_length[kdx]*s_nmk_by_length[my_len-kdx] for kdx in range(kdx_min,kdx_max+1))
        assert counts[my_len] == expected_count
