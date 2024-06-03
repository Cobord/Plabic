"""
Triangulation of convex m-gon
"""

from typing import Tuple,Set,List,Union,cast,Dict,Any
from functools import reduce
import itertools
from math import sin, cos, pi as PI, atan2, sqrt
from .plabic_diagram import PlabicGraph, BiColor
from .framed_2_disks import FramedDiskConfig

Point = Tuple[float,float]

class Triangulation:
    """
    a triangulated convex m-gon
    """

    @staticmethod
    def __makes_triangle(idx_1 : int, idx_2 : int,
                       boxed_off : Set[int],
                       num_points : int) -> Tuple[bool,bool,Union[int,Tuple[int,int]]]:
        """
        if connect points idx_1 and idx_2 does this form a triangle
        after triangles have already been inserted such that the points
        in boxed_off are already behind triangles
        """
        if idx_1==idx_2:
            raise ValueError("Already removed non-diagonals")
        if idx_1>idx_2:
            idx_1,idx_2 = idx_2,idx_1
        whats_between = set(range(idx_1+1,idx_2))
        whats_between_other = set(range(idx_2+1,num_points))
        whats_between_other.update(range(0,idx_1))
        newly_boxed_off_1 = whats_between.difference(boxed_off)
        newly_boxed_off_2 = whats_between_other.difference(boxed_off)
        if len(newly_boxed_off_1)==1:
            is_triangle_1 = True
            box_by_1 = newly_boxed_off_1.pop()
        else:
            is_triangle_1 = False
        if len(newly_boxed_off_2)==1:
            is_triangle_2 = True
            box_by_2 = newly_boxed_off_2.pop()
        else:
            is_triangle_2 = False
        is_triangle = is_triangle_1 or is_triangle_2
        both_triangles = is_triangle_1 and is_triangle_2
        boxed_index : Union[int,Tuple[int,int]] = -1
        if both_triangles:
            #pylint:disable=possibly-used-before-assignment
            boxed_index = (box_by_1,box_by_2)
        elif is_triangle_1:
            boxed_index = box_by_1
        elif is_triangle_2:
            boxed_index = box_by_2
        return is_triangle, both_triangles, boxed_index

    # pylint:disable = too-many-locals, too-many-branches, too-many-statements
    def __init__(self,vertices: List[Point], diagonals : List[Tuple[int,int]]):
        self.num_points = len(vertices)
        self.vertices = vertices
        actually_diagonals = all( (x!=y and 0<=x<self.num_points and
                                   0<=y<self.num_points) for x,y in diagonals)
        if not actually_diagonals:
            raise ValueError("Something was not a diagonal")
        boxed_off : Set[int] = set()
        used_diagonal_idces : Set[int] = set()
        properly_ordered_diagonals : List[Tuple[int,int,Union[int,Tuple[int,int]]]] = []
        for _ in range(len(diagonals)):
            some_diagonal_worked = False
            for idx,cur_diagonal in enumerate(diagonals):
                if idx in used_diagonal_idces:
                    continue
                (from_pt,to_pt) = cur_diagonal
                is_triangle, both_triangles, newly_boxed =\
                    self.__makes_triangle(from_pt, to_pt, boxed_off,self.num_points)
                if both_triangles:
                    some_diagonal_worked = True
                    properly_ordered_diagonals.append((from_pt,to_pt,newly_boxed))
                    break
                if is_triangle:
                    some_diagonal_worked = True
                    properly_ordered_diagonals.append((from_pt,to_pt,newly_boxed))
                    used_diagonal_idces.add(idx)
                    boxed_off.add(cast(int,newly_boxed))
                    break
            if not some_diagonal_worked:
                raise ValueError("No diagonals worked to make a triangle at this stage")
        self.my_diagonal_list = [(x,y) for x,y,_ in properly_ordered_diagonals]
        self.my_quadrilaterals : List[Tuple[int,int,int,int]] = []
        self.my_triangles : List[Tuple[int,int,int]] = []
        for diag_1,diag_2,corner_12 in properly_ordered_diagonals:
            if isinstance(corner_12,tuple):
                self.my_triangles.append((diag_1,corner_12[0],diag_2))
                self.my_triangles.append((diag_1,corner_12[1],diag_2))
                self.my_quadrilaterals.append((diag_1,corner_12[0],diag_2,corner_12[1]))
                continue
            self.my_triangles.append((diag_1,corner_12,diag_2))
            corner_21 = -1
            for diag_3,diag_4,corner_34 in properly_ordered_diagonals:
                if isinstance(corner_34,int):
                    corner_34_cond = corner_34 in [diag_1,diag_2]
                else:
                    corner_34_cond = corner_34[0] in [diag_1,diag_2] or\
                        corner_34[1] in [diag_1,diag_2]
                if corner_34_cond and diag_3 in [diag_1,diag_2]:
                    corner_21 = diag_4
                    break
                if corner_34_cond and diag_4 in [diag_1,diag_2]:
                    corner_21 = diag_3
                    break
            if corner_21<0:
                raise ValueError(f"The other corner for {(diag_1,diag_2,corner_12)} wasn't found")
            self.my_quadrilaterals.append((diag_1,corner_12,diag_2,corner_21))
        aggregated = reduce(lambda acc,x : (acc[0]+x[0],acc[1]+x[1]),self.vertices)
        self.center : Point = (aggregated[0]/self.num_points,aggregated[1]/self.num_points)

    def quad_flip(self,which_diag : Tuple[int,int]) -> Tuple[bool,Tuple[int,int]]:
        """
        flip this diagonal
        returns if such a diagonal existed
        as well as the other diagonal of the quadrilateral
        that is now in the triangulation
        """
        try:
            in_my_diags = self.my_diagonal_list.index(which_diag)
        except ValueError:
            return False, (-1,-1)
        relevant_quad = self.my_quadrilaterals[in_my_diags]
        other_diag = (relevant_quad[1], relevant_quad[3])
        self.my_diagonal_list[in_my_diags] = other_diag
        self.my_quadrilaterals[in_my_diags] = (relevant_quad[1],
                                               relevant_quad[2],
                                               relevant_quad[3],
                                               relevant_quad[0])
        def this_not_a_side(triangle : Tuple[int,int,int]) -> bool:
            """
            check that which_diag is not a side of this triangle
            """
            pt_set = set(triangle)
            pt_set = pt_set.intersection(which_diag)
            return len(pt_set)<2

        self.my_triangles = [cur_triangle
                             for cur_triangle in self.my_triangles
                             if this_not_a_side(cur_triangle)]
        self.my_triangles.append((relevant_quad[1],relevant_quad[2],relevant_quad[3]))
        self.my_triangles.append((relevant_quad[1],relevant_quad[0],relevant_quad[3]))
        return True, (relevant_quad[1],relevant_quad[3])

    def to_frieze(self) -> str:
        """
        The associated SL2 frieze pattern as it is printed out
        """
        num_rows = self.num_points-1
        width_bound = 2*self.num_points
        pattern = [[0 for col in range(width_bound)] for row in range(num_rows)]
        for idx in range(width_bound):
            pattern[0][idx] = 1
            pattern[num_rows-1][idx] = 1
        how_many_triangles = [0]*self.num_points
        for tri in self.my_triangles:
            (ver_a,ver_b,ver_c) = tri
            how_many_triangles[ver_a] += 1
            how_many_triangles[ver_b] += 1
            how_many_triangles[ver_c] += 1
        for idx in range(width_bound):
            pattern[1][idx] = how_many_triangles[idx%self.num_points]
        for row_num in range(2,num_rows-1):
            if row_num%2==0:
                for idx in range(row_num//2-1,width_bound-row_num//2):
                    numerator = pattern[row_num-1][idx]*pattern[row_num-1][idx+1]-1
                    pattern[row_num][idx] = numerator//(pattern[row_num-2][idx])
            else:
                for idx in range(row_num//2,width_bound-row_num//2):
                    numerator = pattern[row_num-1][idx]*pattern[row_num-1][idx-1]-1
                    pattern[row_num][idx] = numerator//(pattern[row_num-2][idx])
        max_num_digits = 1
        for row_num in range(num_rows):
            for (idx,entry) in enumerate(pattern[row_num]):
                max_num_digits = max(max_num_digits,len(str(entry)))
                if entry==0:
                    if idx+self.num_points<width_bound and \
                        pattern[row_num][idx+self.num_points] !=0:
                        pattern[row_num][idx] = pattern[row_num][idx+self.num_points]
                    elif idx>self.num_points and \
                        pattern[row_num][idx-self.num_points] !=0:
                        pattern[row_num][idx] = pattern[row_num][idx-self.num_points]
        spacer = " "*max_num_digits
        questioner = "?"*max_num_digits
        row_strings = (spacer.join([
            questioner if x==0 else str(x)
            for x in pattern[cur_row]])
            for cur_row in range(num_rows))
        return "\n".join((f" {row_string}" if jdx % 2==0 else row_string
                          for jdx,row_string in enumerate(row_strings)))

    def to_plabic(self) -> PlabicGraph:
        """
        a plabic graph whose quiver is the same
        as the quiver associated to triangulation
        """
        my_init_data : Dict[str,Tuple[BiColor,List[str]]] = {}
        extra_node_props : Dict[str,Dict[str,Any]] = {}
        scale_factor = 2
        r2_to_center_all = [0.0]*self.num_points
        for idx in range(self.num_points):
            my_init_data[f"ext{idx}"] = (BiColor.RED,[f"int{idx}"])
            scaled_up_pt = list(self.vertices[idx])
            scaled_up_pt[0] -= self.center[0]
            scaled_up_pt[1] -= self.center[1]
            scaled_up_pt[0] *= scale_factor
            scaled_up_pt[1] *= scale_factor
            r2_to_center = scaled_up_pt[0]**2+scaled_up_pt[1]**2
            r2_to_center_all[idx] = r2_to_center
            scaled_up_pt[0] += self.center[0]
            scaled_up_pt[1] += self.center[1]
            extra_node_props[f"ext{idx}"] = {"position" : (scaled_up_pt[0],scaled_up_pt[1]),
                                             "my_perfect_edge" : -1}
            extra_node_props[f"int{idx}"] = {"position" : self.vertices[idx],
                                             "my_perfect_edge" : 0}
        perfect_matching : Dict[str,str] = {}
        for (triangle_num,(pt_one,pt_two,pt_three)) in enumerate(self.my_triangles):
            loc_one = self.vertices[pt_one]
            loc_two = self.vertices[pt_two]
            loc_three = self.vertices[pt_three]
            averaged_position = (loc_one[0]+loc_two[0]+loc_three[0])/3.0,\
                (loc_one[1]+loc_two[1]+loc_three[1])/3.0
            which_vertices = [f"int{pt_one}",f"int{pt_two}",f"int{pt_three}"]
            which_vertices = clockwise_order(which_vertices,extra_node_props,averaged_position)
            my_matching_pt = f"int{pt_two}"
            where_my_matching_pt = which_vertices.index(my_matching_pt)
            my_init_data[f"triangle{triangle_num}"] =\
                (BiColor.RED,which_vertices)
            extra_node_props[f"triangle{triangle_num}"] = {"position" : averaged_position,
                                                           "my_perfect_edge" : where_my_matching_pt}
            perfect_matching[my_matching_pt] = f"triangle{triangle_num}"
        for idx in range(self.num_points):
            which_vertices = \
                [f"triangle{trinum}" for trinum,tri in enumerate(self.my_triangles)
                 if idx in set(tri)]
            which_vertices.append(f"ext{idx}")
            which_vertices = clockwise_order(which_vertices,extra_node_props,self.vertices[idx])
            my_init_data[f"int{idx}"] = (BiColor.GREEN,which_vertices)
            my_match = perfect_matching.get(f"int{idx}",f"ext{idx}")
            extra_node_props[f"int{idx}"]["my_perfect_edge"] = \
                which_vertices.index(my_match)
            if my_match == f"ext{idx}":
                extra_node_props[f"ext{idx}"]["my_perfect_edge"] = 0
        external_init_orientation = [f"ext{idx}" for idx in range(self.num_points)]
        multi_edge_permutation : Dict[Tuple[str,str],Dict[int,int]] = {}
        internal_bdry_orientations = None
        def get_common_radius(r2_all : List[float]) -> Tuple[bool,float]:
            """
            given bunch of r^2 for all the vertices relative to self.center
            are they all on a circle around that center
            """
            if len(r2_all) < 2:
                raise ValueError("m must be at least 3")
            max_abs_diffs = max(abs(a-b)
                for (a,b) in itertools.combinations(r2_all,2))
            if max_abs_diffs<1e-6*sum(r2_all)/len(r2_all):
                return True,sqrt(r2_all[0])
            return False,0.0
        has_common_radius, common_radius = get_common_radius(r2_to_center_all)
        if has_common_radius:
            my_circles_config = FramedDiskConfig([],(self.center,common_radius,0.0))
        else:
            my_circles_config = None
        return PlabicGraph(my_init_data,
                 external_init_orientation,
                 multi_edge_permutation,
                 internal_bdry_orientations,
                 extra_node_props,
                 my_circles_config)

def clockwise_order(node_names : List[str],
                    prop_dict : Dict[str,Dict[str,Any]],
                    center : Point) -> List[str]:
    """
    order the node_names clockwise with the usual cut
    based on their positions in prop_dict
    relative to the center
    """
    node_positions : List[Point] = [prop_dict[cur_pt]["position"] for cur_pt in node_names]
    node_vectors = [(x-center[0],y-center[0]) for x,y in node_positions]
    node_name_angle = list(zip(node_names,[atan2(y,x) for (x,y) in node_vectors]))
    node_name_angle.sort(key=lambda z: -z[1])
    return [node_name for (node_name,_) in node_name_angle]

if __name__ == '__main__':
    N = 8
    my_diagonals = [(4,7),(2,8),(2,4),(2,7),(4,6)]
    RADIUS = 2
    t = Triangulation([ (RADIUS*cos(2*PI*which/N),RADIUS*sin(2*PI*which/N))
                       for which in range(N)], [(x-1,y-1) for x,y in my_diagonals])
    p = t.to_plabic()
    p.draw()
    for tri_num in range(6):
        cur_triangle_name = f"triangle{tri_num}"
        assert p.my_kn_type(cur_triangle_name) == (1,8)
    for int_num,exp_contrib in enumerate([2,5,2,5,2,3,4,3]):
        cur_int_name = f"int{int_num}"
        assert p.my_kn_type(cur_int_name) == (2-exp_contrib,8)
    assert p.my_kn_type() == (2,8),f"{p.my_kn_type()}"
    change, new_diag = t.quad_flip((1,3))
    print(f"Triangles are : {t.my_triangles}")
    p = t.to_plabic()
    assert p.my_kn_type() == (2,8),f"{p.my_kn_type()}"
    print(f"Perfect matching : {p.my_perfect_matching}")
    p.draw()

    N = 5
    my_diagonals = [(1,3),(1,4)]
    RADIUS = 2
    t = Triangulation([ (RADIUS*cos(2*PI*which/N),RADIUS*sin(2*PI*which/N))
                       for which in range(N)], [(x-1,y-1) for x,y in my_diagonals])
    print(t.to_frieze())

    N = 8
    my_diagonals = [(8,6),(8,5),(5,3),(1,3),(8,3)]
    RADIUS = 4
    t = Triangulation([ (RADIUS*cos(2*PI*which/N),RADIUS*sin(2*PI*which/N))
                       for which in range(N)], [(x-1,y-1) for x,y in my_diagonals])
    print(t.to_frieze())
