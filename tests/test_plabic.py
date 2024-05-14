"""
test Plabic graphs
"""
#pylint:disable=too-many-locals
from typing import List,Dict
from plabic import PlabicGraph,BiColor,PlabicGraphBuilder

def test_plabic() -> None:
    """
    example from figure 7.1 of https://arxiv.org/pdf/2106.02160.pdf
    black there -> red here
    white there -> green here
    the names of internal vertices are numbered top to bottom
    left to right and noting the one attached to external 6
    is slightly below the 3rd and 4th (though that is hard to see in the picture)
    """
    external_orientation = ["ext1", "ext2", "ext3", "ext4", "ext5", "ext6"]
    _internal_vertices = ["int1", "int2", "int3",
                         "int4", "int5", "int6", "int7", "int8"]
    my_data = {}
    my_data["ext1"] = (BiColor.RED, ["int1"])
    my_data["ext2"] = (BiColor.RED, ["int2"])
    my_data["ext3"] = (BiColor.RED, ["int7"])
    my_data["ext4"] = (BiColor.RED, ["int8"])
    my_data["ext5"] = (BiColor.RED, ["int6"])
    my_data["ext6"] = (BiColor.RED, ["int5"])

    my_data["int1"] = (BiColor.GREEN, ["ext1", "int2", "int3"])
    my_data["int2"] = (BiColor.RED, ["int1", "ext2", "int4"])
    my_data["int3"] = (BiColor.RED, ["int1", "int4", "int6"])
    my_data["int4"] = (BiColor.GREEN, ["int2", "int7", "int3"])
    my_data["int5"] = (BiColor.GREEN, ["ext6"])
    my_data["int6"] = (BiColor.GREEN, ["int3", "int8", "ext5"])
    my_data["int7"] = (BiColor.GREEN, ["int4", "ext3", "int8"])
    my_data["int8"] = (BiColor.RED, ["int6", "int7", "ext4"])
    example_plabic = PlabicGraph(my_data, external_orientation, {})
    assert example_plabic.my_extra_props.issubset([])
    expected_one_steps = ["int1", "int2", "int7", "int8", "int6", "int5"]
    expected_two_steps = ["int2", "int1", "int8", "int7", "int3", "ext6"]
    for cur_bdry_temp, exp_one_step, exp_two_step in\
            zip(external_orientation, expected_one_steps, expected_two_steps):
        one_step = example_plabic.follow_rules_of_road(cur_bdry_temp, None, None)
        assert one_step is not None
        cur_one_step_idx, cur_one_step = one_step
        assert cur_one_step == exp_one_step
        two_step = example_plabic.follow_rules_of_road(
            cur_one_step, cur_bdry_temp, cur_one_step_idx)
        assert two_step is not None
        assert two_step[1] == exp_two_step
    expected_perm = {"ext1":"ext3", "ext2":"ext4",
                         "ext3":"ext5","ext4":"ext1",
                         "ext5":"ext2","ext6":"ext6"}
    expected_decorated = {"ext6":BiColor.GREEN}
    correct_decorated_permutation(example_plabic,external_orientation,
                                  expected_perm,expected_decorated)
    changed, explanation = example_plabic.square_move(("int1", "int2", "int4", "int3"))
    assert changed and explanation == "Success"
    correct_decorated_permutation(example_plabic,external_orientation,
                                  expected_perm,expected_decorated)
    changed, explanation = example_plabic.square_move(("int1", "int2", "int4", "int3"))
    assert changed and explanation == "Success"
    correct_decorated_permutation(example_plabic,external_orientation,
                                  expected_perm,expected_decorated)
    changed, explanation = example_plabic.flip_move("int4", "int7")
    assert changed and explanation == "Success"
    correct_decorated_permutation(example_plabic,external_orientation,
                                  expected_perm,expected_decorated)
    simplifies = example_plabic.greedy_shrink()
    assert simplifies
    simplifies = example_plabic.greedy_shrink()
    assert not simplifies
    correct_decorated_permutation(example_plabic,external_orientation,
                                  expected_perm,expected_decorated)
    did_scale_up = example_plabic.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert not did_scale_up

# pylint:disable = too-many-statements
def test_builder():
    """
    same example as before but with an internal circle
    which has 1 vertex which connects to int5
    """
    my_builder = PlabicGraphBuilder()
    my_builder.set_num_external(6)
    my_builder.add_external_bdry_vertex("ext1",0,"int1")
    my_builder.add_external_bdry_vertex("ext2",1,"int2")
    my_builder.add_external_bdry_vertex("ext3",2,"int7")
    my_builder.add_external_bdry_vertex("ext4",3,"int8")
    my_builder.add_external_bdry_vertex("ext5",4,"int6")
    my_builder.add_external_bdry_vertex("ext6",5,"int5")
    my_builder.set_internal_circles_nums([1])
    my_builder.add_internal_bdry_vertex("bdry?",0,0,"int5")
    my_builder.add_internal_vertex("int1",BiColor.GREEN, ["ext1", "int2", "int3"])
    my_builder.add_internal_vertex("int2",BiColor.RED, ["int1", "ext2", "int4"])
    my_builder.add_internal_vertex("int3",BiColor.RED, ["int1", "int4", "int6"])
    my_builder.add_internal_vertex("int4",BiColor.GREEN, ["int2", "int7", "int3"])
    my_builder.add_internal_vertex("int5",BiColor.GREEN, ["ext6","bdry?"])
    my_builder.add_internal_vertex("int6",BiColor.GREEN, ["int3", "int8", "ext5"])
    my_builder.add_internal_vertex("int7",BiColor.GREEN, ["int4", "ext3", "int8"])
    my_builder.add_internal_vertex("int8",BiColor.RED, ["int6", "int7", "ext4"])
    my_builder.rename("bdry?","bdry1")
    example_plabic = my_builder.build()
    assert example_plabic.my_extra_props.issubset([])
    exp_external_orientation = ["ext1", "ext2", "ext3", "ext4", "ext5", "ext6"]
    assert example_plabic.my_external_nodes == exp_external_orientation
    expected_one_steps = ["int1", "int2", "int7", "int8", "int6", "int5"]
    expected_two_steps = ["int2", "int1", "int8", "int7", "int3", "bdry1"]
    for cur_bdry_temp, exp_one_step, exp_two_step in\
            zip(exp_external_orientation, expected_one_steps, expected_two_steps):
        one_step = example_plabic.follow_rules_of_road(cur_bdry_temp, None, None)
        assert one_step is not None
        cur_one_step_idx, cur_one_step = one_step
        assert cur_one_step == exp_one_step
        two_step = example_plabic.follow_rules_of_road(
            cur_one_step, cur_bdry_temp, cur_one_step_idx)
        assert two_step is not None
        assert two_step[1] == exp_two_step
    expected_perm = {"ext1":"ext3", "ext2":"ext4",
                         "ext3":"ext5","ext4":"ext1",
                         "ext5":"ext2","ext6":"bdry1","bdry1":"ext6"}
    expected_decorated = {}
    correct_decorated_permutation(example_plabic,exp_external_orientation,
                                  expected_perm,expected_decorated)
    changed, explanation = example_plabic.square_move(("int1", "int2", "int4", "int3"))
    assert changed and explanation == "Success"
    correct_decorated_permutation(example_plabic,exp_external_orientation,
                                  expected_perm,expected_decorated)
    changed, explanation = example_plabic.square_move(("int1", "int2", "int4", "int3"))
    assert changed and explanation == "Success"
    correct_decorated_permutation(example_plabic,exp_external_orientation,
                                  expected_perm,expected_decorated)
    changed, explanation = example_plabic.flip_move("int4", "int7")
    assert changed and explanation == "Success"
    correct_decorated_permutation(example_plabic,exp_external_orientation,
                                  expected_perm,expected_decorated)
    changed, explanation = example_plabic.remove_bivalent("int5")
    assert not changed and explanation == "Cannot cause two boundary nodes to connect directly"
    simplifies = example_plabic.greedy_shrink()
    assert simplifies
    simplifies = example_plabic.greedy_shrink()
    assert not simplifies
    correct_decorated_permutation(example_plabic,exp_external_orientation,
                                  expected_perm,expected_decorated)
    did_scale_up = example_plabic.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert not did_scale_up

def correct_decorated_permutation(example_plabic : PlabicGraph,
                                        external_orientation : List[str],
                                        expected_perm : Dict[str,str],
                                        expected_decorated : Dict[str,BiColor]):
    """
    assert that the decorated permutation is as expected
    """
    for cur_bdry_temp in external_orientation:
        turn_color, cur_tgt = example_plabic.bdry_to_bdry(cur_bdry_temp)
        assert cur_tgt == expected_perm[cur_bdry_temp]
        assert turn_color == expected_decorated.get(cur_bdry_temp,None)
