"""
test for (bounded) affine permutations
"""
import itertools
import random
from math import factorial
from typing import Dict,Set
from plabic import AffinePermutation, BoundedAffinePermutation, BruhatInterval, BiColor

def test_coxeter_muls():
    """
    test coxeter
    """
    trial_n = 5
    cox_gens = {}
    for idx in range(trial_n):
        s_idx = AffinePermutation(n_val=trial_n,coxeter_word=[idx])
        assert s_idx.to_reduced_word() == [idx]
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
        assert combined_exp.to_reduced_word() == [idx,jdx]
        assert combined_exp.my_coxeter_word == [idx,jdx]
        assert combined_obs.my_coxeter_word == [idx,jdx]
        assert combined_exp == combined_obs
    for (idx,jdx) in itertools.permutations(range(trial_n),2):
        s_idx = cox_gens[idx]
        s_jdx = cox_gens[jdx]
        if trial_n-1>abs(idx-jdx)>1:
            combined_exp = AffinePermutation(n_val=trial_n,coxeter_word=[jdx])
            expected_reduced_word_exp = [jdx]
            shrinks = True
        else:
            combined_exp = AffinePermutation(n_val=trial_n,coxeter_word=[jdx,idx,jdx])
            expected_reduced_word_exp = [jdx,idx,jdx]
            shrinks = False
        combined_obs = s_idx*s_jdx*s_idx
        assert combined_obs.my_coxeter_word == [idx,jdx,idx]
        assert combined_exp.my_coxeter_word == expected_reduced_word_exp
        assert combined_exp == combined_obs
        if shrinks:
            assert combined_obs.my_coxeter_word == combined_exp.my_coxeter_word
        else:
            assert combined_obs.my_coxeter_word == [idx,jdx,idx]
    for (idx,jdx,kdx) in itertools.permutations(range(trial_n),3):
        s_idx = cox_gens[idx]
        s_jdx = cox_gens[jdx]
        s_kdx = cox_gens[kdx]
        combined_obs = AffinePermutation(n_val=trial_n,coxeter_word=[idx,jdx,kdx])
        combined_exp = s_idx*s_jdx*s_kdx
        assert combined_exp.to_reduced_word() == [idx,jdx,kdx]
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

def test_reduced():
    """
    test that to_reduced_word produces an equivalent element
    """
    my_n = 9
    trials_num = 30
    word_lengths = [0,1,2,3,5,10]
    for _ in range(trials_num):
        random_word = random.choices(range(my_n), k = random.choice(word_lengths))
        word_initial = AffinePermutation(coxeter_word=random_word,n_val=my_n)
        assert word_initial.my_coxeter_word == random_word
        reduced_random_word = word_initial.to_reduced_word()
        word_reconstructed = AffinePermutation(coxeter_word=reduced_random_word,n_val=my_n)
        assert len(reduced_random_word)<=len(random_word)
        if len(reduced_random_word)<len(random_word):
            assert word_initial.my_coxeter_word == reduced_random_word
        else:
            assert word_initial.my_coxeter_word == random_word
        assert word_initial == word_reconstructed
        assert word_reconstructed.my_coxeter_word == reduced_random_word
        if len(reduced_random_word)<len(random_word):
            assert word_initial.my_coxeter_word == reduced_random_word
        else:
            assert word_initial.my_coxeter_word == random_word

def test_jumpers():
    """
    ij_jumpers which is a helper for the Bruhat order
    """
    test_perm = AffinePermutation(some_vals = {1:2,2:0,3:4}, n_val = 3)
    assert test_perm.ij_jumpers(3,1) == {3,1,0}

def test_qkn():
    """
    Q_kn
    """
    my_k = 2
    my_n = 7
    max_length = 5
    for _cur_interval in BruhatInterval.all_Qkn(my_k,my_n,max_length):
        pass

def test_bruhat():
    """
    see if using the sets for bruhat_leq
    """
    perm_1 = AffinePermutation(some_vals = {1:1,2:2,3:3}, n_val = 3)
    perm_2 = AffinePermutation(some_vals = {1:2,2:1,3:3}, n_val = 3)
    perm_3 = AffinePermutation(some_vals = {1:2,2:0,3:4}, n_val = 3)
    known_leq_dict : Dict[AffinePermutation,Set[AffinePermutation]] = {} #type:ignore
    known_geq_dict : Dict[AffinePermutation,Set[AffinePermutation]] = {} #type:ignore
    known_leq_dict[perm_1] = set()
    known_leq_dict[perm_2] = set()
    known_leq_dict[perm_3] = set()
    known_geq_dict[perm_1] = set()
    known_geq_dict[perm_2] = set()
    known_geq_dict[perm_3] = set()
    assert len(known_leq_dict[perm_2]) == 0
    assert len(known_geq_dict[perm_1]) == 0
    assert perm_1.bruhat_leq(perm_2,
                             known_geq_dict.get(perm_1,None),
                             known_leq_dict.get(perm_2,None))
    assert len(known_leq_dict[perm_2]) == 1
    assert len(known_geq_dict[perm_1]) == 1
    assert perm_1 in known_leq_dict[perm_2]
    assert perm_2 in known_geq_dict[perm_1]
    assert perm_1.bruhat_leq(perm_2,
                             known_geq_dict.get(perm_1,None),
                             known_leq_dict.get(perm_2,None))
    assert len(known_leq_dict[perm_2]) == 1
    assert len(known_geq_dict[perm_1]) == 1
    assert perm_1 in known_leq_dict[perm_2]
    assert perm_2 in known_geq_dict[perm_1]
    assert perm_2.bruhat_leq(perm_3,
                             known_geq_dict.get(perm_2,None),
                             known_leq_dict.get(perm_3,None))
    assert len(known_leq_dict[perm_3]) == 1
    assert len(known_geq_dict[perm_2]) == 1
    assert perm_2 in known_leq_dict[perm_3]
    assert perm_3 in known_geq_dict[perm_2]
    assert perm_1.bruhat_leq(perm_3,
                             known_geq_dict.get(perm_1,None),
                             known_leq_dict.get(perm_3,None))
    assert len(known_leq_dict[perm_3]) == 2
    assert len(known_geq_dict[perm_1]) == 2
    assert perm_1 in known_leq_dict[perm_3]
    assert perm_2 in known_leq_dict[perm_3]
    assert perm_3 in known_geq_dict[perm_1]
    assert perm_2 in known_geq_dict[perm_1]

def test_plabic():
    """
    plabic graphs constructed using bridge decomposition
    """
    my_n = 4
    for cur_ba in BoundedAffinePermutation.all_bounded_affine_perms(my_n):
        cur_plabic = cur_ba.to_plabic()
        assert "position" in cur_plabic.my_extra_props
        for idx in range(1,my_n+1):
            fix_pt_info,idx_goes_to = cur_plabic.bdry_to_bdry(f"ext{idx}")
            expected_go_to_num = cur_ba[idx] % my_n if cur_ba[idx] % my_n != 0 else my_n
            if fix_pt_info is not None:
                if fix_pt_info == BiColor.GREEN:
                    assert cur_ba[idx]==idx+my_n
                else:
                    assert cur_ba[idx]==idx
            assert idx_goes_to == f"ext{expected_go_to_num}"
