"""
a framed little 2-disk configuration
"""
from __future__ import annotations
from typing import Callable,Tuple,List,Optional,cast
import numpy as np
import matplotlib.pyplot as plt

Point = Tuple[float,float]
Radius = float
Angle = float

class FramedDiskConfig:
    """
    a framed little 2-disk configuration
    """

    def __init__(self,
                 internal_circles : List[Tuple[Point,Radius,Angle]],
                 outer_circle : Optional[Tuple[Point,Radius,Angle]] = None):
        if outer_circle is None:
            outer_circle = ((0.0,0.0),1.0,0.0)
        self.outer_circle = outer_circle
        self.internal_circles = internal_circles.copy()

    # pylint:disable=too-many-locals
    def do_transformation(self,which_disk : int,
                           other : FramedDiskConfig) -> Callable[[Point],Point]:
        """
        the transformation necessary to take which_disk internal of self
        to the outer_circle of other
        and does that on self
        """
        if which_disk < 0:
            raise ValueError("Which disk is a natural number")
        if len(self.internal_circles) <= which_disk:
            raise ValueError("Not enough internal disks")
        target_center, target_radius, target_angle = other.outer_circle
        cur_center, cur_radius, cur_angle = self.internal_circles[which_disk]
        as_np_cur_center = np.array(cur_center)
        scale_factor = target_radius/cur_radius
        angle_diff = target_angle-cur_angle
        cos_angle_diff, sin_angle_diff = \
            np.cos(angle_diff), np.sin(angle_diff)
        rot_matrix = np.array(((cos_angle_diff, -sin_angle_diff),
                                (sin_angle_diff, cos_angle_diff)))
        as_np_target_center = np.array(target_center)
        def point_transformation(my_point : Point) -> Point:
            as_np_mine = np.array(my_point)
            centered_pt = (as_np_mine - as_np_cur_center)*scale_factor
            centered_rotated_pt = np.matmul(rot_matrix,centered_pt)
            final_pt = centered_rotated_pt + as_np_target_center
            return cast(Point,tuple(final_pt.tolist()))
        def circle_transformation(framed_circle : Tuple[Point,Radius,Angle])\
            -> Tuple[Point,Radius,Angle]:
            (old_center,old_radius,old_angle) = framed_circle
            new_center = point_transformation(old_center)
            assert isinstance(new_center,tuple) and len(new_center)==2 and\
                isinstance(new_center[0],float) and isinstance(new_center[1],float)
            return (new_center,
                    old_radius*scale_factor,
                    old_angle+angle_diff)
        self.outer_circle = circle_transformation(self.outer_circle)
        _ = self.internal_circles.pop(which_disk)
        self.internal_circles = [circle_transformation(z) for z in self.internal_circles]
        self.internal_circles.extend(other.internal_circles)
        return point_transformation

    #pylint:disable=too-many-arguments
    def draw(self,on_top_of_existing : bool = False,
             show_as_well : bool = True,
             draw_framing : bool = True,
             external_circle_color : str="black",
             internal_circles_color : str="black"):
        """
        draw circles
        """
        if on_top_of_existing:
            fig = plt.gcf()
            my_ax = fig.gca()
        else:
            fig, my_ax = plt.subplots()
        (center_x,center_y),radius,_ = self.outer_circle
        if len(self.internal_circles)==0 or not draw_framing:
            min_radius = radius
        else:
            min_radius = min(internal_radius for _,internal_radius,_ in self.internal_circles)
        dot_radius = min_radius / 10
        plt.ylim(center_y-radius-2*dot_radius,center_y+radius+2*dot_radius)
        plt.xlim(center_x-radius-2*dot_radius,center_x+radius+2*dot_radius)

        def draw_a_circle(framed_circle : Tuple[Point,Radius,Angle],color : str):
            (center_x,center_y),radius,angle = framed_circle
            circle = plt.Circle((center_x, center_y), radius, color=color,fill=False)
            my_ax.add_patch(circle)
            if draw_framing:
                cos_angle, sin_angle = \
                    np.cos(angle), np.sin(angle)
                rot_matrix = np.array(((cos_angle, -sin_angle),
                                        (sin_angle, cos_angle)))
                marked_point = np.array([center_x,center_y])+\
                    radius*np.matmul(rot_matrix,np.array((1.0,0.0)))
                marked_point_list = marked_point.tolist()
                xmarked,ymarked = marked_point_list[0], marked_point_list[1]
                circle2 = plt.Circle((xmarked, ymarked), dot_radius, color=color,fill=True)
                my_ax.add_patch(circle2)

        draw_a_circle(self.outer_circle,external_circle_color)
        for cur_internal in self.internal_circles:
            draw_a_circle(cur_internal,internal_circles_color)

        plt.draw()
        if not on_top_of_existing and show_as_well:
            plt.show()

if __name__ == '__main__':
    e = FramedDiskConfig([])
    e.draw()

    e1 = FramedDiskConfig([((0,.5),.2,0)])
    e1.draw()

    e2 = FramedDiskConfig([((0,.5),.2,0),((-.5,0),.3,3.14)])
    e2.draw()

    _ = e1.do_transformation(0,e2)
    e1.draw(draw_framing=False)

    _ = e1.do_transformation(0,e2)
    e1.draw(draw_framing=False)
