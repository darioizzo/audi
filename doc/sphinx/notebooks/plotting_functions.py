from pyaudi import gdual_double as gdual, taylor_model, int_d, generate_combinations
import numpy as np
import seaborn as sns
import matplotlib.colors as col
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from scipy.integrate import solve_ivp
import copy
import matplotlib
from collections import Counter
import trimesh


def plot_taylor_model_1d(
    tm: taylor_model,
    ax=None,
    n_of_points=100,
    buffer=0.1,
    color="b",
):
    """
    Plot a 1-D Taylor model.

    Args:
        tm (taylor_model): Taylor model object.
        ax (matplotlib.axes.Axes): The axes to plot on.
        n_of_points (int): The number of points to plot.
        buffer (float): The buffer to use for the plot.
        color (str): The color of the plot.

    Returns:
        matplotlib.axes.Axes: The axes with the plotted Taylor polynomial.
    """

    # Get plotting points
    var, dom = next(iter(tm.domain.items()))
    plotting_points = np.linspace(dom.lower, dom.upper, n_of_points)

    # Evaluate shifted polynomial
    shifted_func_values = np.zeros_like(plotting_points)
    for i, pt in enumerate(plotting_points):
        shifted_func_values[i] = tm.tpol.evaluate({"d" + var: (pt - tm.exp_point[var])})

    # Get upper and lower bound values
    lower_bound_func_values = shifted_func_values + tm.rem_bound.lower
    upper_bound_func_values = shifted_func_values + tm.rem_bound.upper

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    ax.plot(plotting_points, lower_bound_func_values, c="k", linewidth=0.5)
    ax.plot(plotting_points, upper_bound_func_values, c="k", linewidth=0.5)
    ax.fill_between(
        plotting_points, lower_bound_func_values, upper_bound_func_values, color=color
    )

    ax.set_xlim([dom.lower - buffer, dom.upper + buffer])

    return ax


def plot_taylor_model_2d(
    tm: taylor_model,
    ax=None,
    n_of_points=100,
    buffer=0.1,
    color="b",
    plot_type="surface",
):
    """
    Plot a 2-D Taylor model.

    Args:
        tm (taylor_model): Taylor model object.
        ax (matplotlib.axes.Axes): The axes to plot on.
        n_of_points (int): The number of points to plot.
        buffer (float): The buffer to use for the plot.
        color (str): The color of the plot.
        plot_type (str): The type of plot. Options are "surface", "contour", or "wireframe".

    Returns:
        matplotlib.axes.Axes: The axes with the plotted Taylor polynomial.
    """

    # Get plotting points
    plotting_points = {}
    for var, dom in tm.domain.items():
        plotting_points[var] = np.linspace(dom.lower, dom.upper, n_of_points)

    # Evaluate shifted polynomial
    domain_iter = iter(tm.domain.items())
    sym1, domain_1 = next(domain_iter)
    sym2, domain_2 = next(domain_iter)
    X, Y = np.meshgrid(
        np.linspace(domain_1.lower, domain_1.upper, n_of_points),
        np.linspace(domain_2.lower, domain_2.upper, n_of_points),
    )

    shifted_func_values = np.array(
        [
            [
                tm.tpol.evaluate(
                    {
                        "d" + sym1: (x - tm.exp_point[sym1]),
                        "d" + sym2: (y - tm.exp_point[sym2]),
                    }
                )
                for x, y in zip(row_x, row_y)
            ]
            for row_x, row_y in zip(X, Y)
        ]
    )

    # Get upper and lower bound values
    lower_bound_func_values = shifted_func_values + tm.rem_bound.lower
    upper_bound_func_values = shifted_func_values + tm.rem_bound.upper

    if ax is None:
        fig = plt.figure(figsize=(7.5, 5))
        ax = fig.add_subplot(111)
    if plot_type == "surface":
        ax.plot_surface(X, Y, lower_bound_func_values, alpha=0.4, color=color)
        ax.plot_surface(X, Y, upper_bound_func_values, alpha=0.4, color=color)
    elif plot_type == "contour":
        ax.contourf(
            X, Y, lower_bound_func_values, levels=n_of_points, cmap="gray", alpha=0.5
        )
        ax.contourf(
            X, Y, upper_bound_func_values, levels=n_of_points, cmap="gray", alpha=0.5
        )
    elif plot_type == "wireframe":
        ax.plot_wireframe(X, Y, lower_bound_func_values, color="k", alpha=1.0)
        ax.plot_wireframe(X, Y, upper_bound_func_values, color="k", alpha=1.0)

    ax.set_xlim([domain_1.lower - buffer, domain_1.upper + buffer])
    ax.set_ylim([domain_2.lower - buffer, domain_2.upper + buffer])
    ax.set_zlim(
        [
            min(lower_bound_func_values.flatten()) - buffer,
            max(upper_bound_func_values.flatten()) + buffer,
        ]
    )

    return ax


def plot_func(ax, func, domain: dict[str, int_d], plot_type="surface", n_of_points=100):
    """
    Plot a function in 1D or 2D. This function is to visualize the actual function that is being approximated by the Taylor model.

    Args:
        ax: The axes to plot on.
        func: The function to plot.
        domain: The domain of the function.
        plot_type: The type of plot (surface, contour, wireframe).
        n_of_points: The number of points to use for plotting.

    Returns:
        The updated axes with the plotted function.
    """
    if isinstance(domain, dict) and len(domain) == 1:
        return plot_func_1d(ax, func, domain)
    elif isinstance(domain, dict) and len(domain) == 2:
        return plot_func_2d(
            ax, func, domain, plot_type=plot_type, n_of_points=n_of_points
        )
    else:
        raise RuntimeError("Plotting an n-dimensional surface is not implemented.")


def plot_func_1d(ax, func, domain: dict[str, int_d]):
    """
    Plot a function in 1D.

    Args:
        ax: The axes to plot on.
        func: The function to plot.
        domain: The domain of the function.

    Returns:
        The updated axes with the plotted function.
    """
    domain = next(iter(domain.values()))
    x_vals = np.linspace(domain.lower, domain.upper, 100)
    y_vals = [func(i) for i in x_vals]
    ax.plot(x_vals, y_vals, c="k")
    return ax


def plot_func_2d(ax, func, domain: dict[int_d], plot_type="surface", n_of_points=100):
    """
    Plot a function in 2D.

    Args:
        ax: The axes to plot on.
        func: The function to plot.
        domain: The domain of the function.
        plot_type: The type of plot (surface, contour, wireframe).
        n_of_points: The number of points to use for plotting.

    Returns:
        The updated axes with the plotted function.
    """
    x_vals = np.linspace(
        float(domain["x"].lower), float(domain["x"].upper), n_of_points
    )
    y_vals = np.linspace(
        float(domain["y"].lower), float(domain["y"].upper), n_of_points
    )
    X, Y = np.meshgrid(x_vals, y_vals)
    Z = np.array(
        [[func(x, y) for x, y in zip(row_x, row_y)] for row_x, row_y in zip(X, Y)]
    )

    if plot_type == "surface":
        ax.plot_surface(X, Y, Z, color="k", alpha=1.0)
    elif plot_type == "contour":
        ax.contourf(X, Y, Z, levels=n_of_points, cmap="gray", alpha=0.5)
    elif plot_type == "wireframe":
        ax.plot_wireframe(X, Y, Z, color="k", alpha=1.0)

    return ax


def add_to_legend(ax, color, label, is_patch=True, loc="best"):
    """
    Add a label to the legend of a plot.

    Args:
        ax: The axes to plot on.
        color: The color of the label.
        label: The label to add.
        is_patch: Whether the label is a patch or not.
        loc: The location of the legend.

    Returns:
        The updated axes with the added label in the legend.
    """
    # Get existing legend
    legend = ax.get_legend()
    if legend is not None:
        existing_handles, existing_labels = legend.legend_handles, [
            t.get_text() for t in legend.texts
        ]
    else:
        existing_handles, existing_labels = [], []

    # Add new legend entry
    if label not in existing_labels:  # Avoid duplicates
        if is_patch:
            existing_handles.append(Patch(color=color, label=label))
        else:
            existing_handles.append(Line2D([0], [0], color=color, label=label, lw=1))
        existing_labels.append(label)

    # Update legend
    ax.legend(handles=existing_handles, loc=loc)

    return ax


def plot_orders(
    func,
    rem_bound,
    exp_point,
    domain,
    ax=None,
    orders=[2, 3],
    buffer=0.1,
    plot_type="surface",
    n_of_points=100,
):
    """
    Plot the Taylor polynomial of a function with given orders. Orders should be in ascending order.

    Args:
        ax: The axes to plot on.
        func: The function to plot.
        rem_bound: The remainder bound.
        exp_point: The expansion point.
        domain: The domain of the function.
        orders: The orders of the Taylor polynomial.
        buffer: The buffer for the plot.
        plot_type: The type of plot (surface, contour, wireframe).
        n_of_points: The number of points to use for plotting.

    Returns:
        The updated axes with the plotted Taylor polynomial and orders.
    """

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    # Plot func
    ax = plot_func(ax, func, domain, plot_type=plot_type, n_of_points=n_of_points)
    # fx_label = "f(x, y)" if isinstance(ax, Axes3D) else "f(x)"
    fx_label = "f(x)"
    ax = add_to_legend(ax, "#000000", fx_label)

    # Plot orders of Taylor model
    for it, order in enumerate(orders):
        colors = sns.color_palette("Blues", as_cmap=True)
        rgb_color = col.rgb2hex(colors(100 * (it + 1)))

        # For 2-D
        if (
            isinstance(exp_point, dict)
            and isinstance(domain, dict)
            and len(exp_point) == 2
            and len(domain) == 2
        ):

            domain_iter = iter(domain.items())
            sym1, domain_1 = next(domain_iter)
            sym2, domain_2 = next(domain_iter)
            x = gdual(
                exp_point[sym1], sym1, order
            )  # Centered at 0.0, called "x1", to order 5
            T_x = taylor_model(
                x,
                rem_bound,
                {k: v for k, v in exp_point.items() if k in [sym1]},
                {k: v for k, v in domain.items() if k in [sym1]},
            )
            y = gdual(
                exp_point[sym2], sym2, order
            )  # Centered at 0.0, called "x1", to order 5
            T_y = taylor_model(
                y,
                rem_bound,
                {k: v for k, v in exp_point.items() if k in [sym2]},
                {k: v for k, v in domain.items() if k in [sym2]},
            )
            f_T = func(T_x, T_y)
            ax = plot_taylor_model_2d(
                f_T,
                ax,
                color=rgb_color,
                buffer=buffer,
                n_of_points=n_of_points,
                plot_type=plot_type,
            )
        # For 1-D
        elif (
            isinstance(exp_point, dict)
            and isinstance(domain, dict)
            and len(exp_point) == 1
            and len(domain) == 1
        ):
            x = gdual(
                exp_point["x"], "x", order
            )  # Centered at 0.0, called "x1", to order 5
            T_x = taylor_model(x, rem_bound, exp_point, domain)
            f_T = func(T_x)
            ax = plot_taylor_model_1d(
                f_T,
                ax,
                color=rgb_color,
                buffer=buffer,
                n_of_points=n_of_points,
            )
        else:
            raise RuntimeError(
                "The input types are not correct. Review the dimensions of your inputs."
            )

        ax = add_to_legend(ax, rgb_color, f"Order: {order}")

    return ax


##############################
### For enclosure plotting ###
##############################


def get_enclosure_line_coords_ndim(
    tm,
    symb,
    symb_dict,
    resolution=100,
):
    """
    Get the coordinates of the enclosure line for a given variable index.

    Args:
        tm (taylor_model): The Taylor model.
        symb (str): The variable.
        symb_dict (dict): Dict containing variables assigned an '1' or '0'
        depending on whether to take upper or lower bound.
        resolution (int): The number of points to use for plotting.

    Returns:
        np.ndarray: The coordinates of the enclosure line.
    """

    # Get variable domains
    var_symb_d = tm.domain[symb]
    var_symb_ds = np.abs(var_symb_d.upper - var_symb_d.lower)

    # Get domain values of constant variables
    const_symb_values = {}
    for const_symb, ind in symb_dict.items():
        if ind == 1:
            const_symb_values[const_symb] = tm.domain[const_symb].upper
        elif ind == 0:
            const_symb_values[const_symb] = tm.domain[const_symb].lower
        else:
            raise RuntimeError(
                f"Invalid entry of {const_minmax_indices}. Only 0 and 1 allowed."
            )

    # Assign dict of const values to use for evaluation
    evaluate_dict = {
        "d" + temp_symb: (const_symb_values[temp_symb] - tm.exp_point[temp_symb])
        for temp_symb in symb_dict
    }

    # Create list of coordinates of variable and evaluate polynomial
    var_symb_vals = np.linspace(-var_symb_ds / 2, var_symb_ds / 2, resolution)
    symb_vals = []
    for val in var_symb_vals:
        evaluate_dict["d" + symb] = val
        symb_vals.append(tm.tpol.evaluate(evaluate_dict))

    return symb_vals


def get_solution_enclosure_ndim(
    tm: taylor_model,
    resolution=100,
):
    """
    Get the solution enclosure between two Taylor models in n dimensions.
    Can have n variables as long as it can be mapped to two coordinates in the end

    Args:
        tm (taylor_model): The first Taylor model.
        resolution (int): The number of points to use for plotting.

    Returns:
        tuple: A tuple containing the center point and the list of coordinates.
    """

    cp1 = tm.tpol.evaluate(
        {"d" + tm.tpol.symbol_set[i]: 0.0 for i in range(tm.ndim)}
    )  # Deviation formulation is formulated around x = 0.0

    hypercube_points = generate_combinations(
        [1 for _ in range(tm.ndim - 1)], cap_sum_indices=False
    )

    coords_list = []
    # Loop over all variables as the 'independent' variable'
    for symb in tm.tpol.symbol_set:

        other_symbs = copy.deepcopy(tm.tpol.symbol_set)
        other_symbs.remove(symb)
        list_of_other_symb_dicts = []
        for hyper_point in hypercube_points:
            list_of_other_symb_dicts.append(
                {other_symbs[ind]: onoff for ind, onoff in enumerate(hyper_point)}
            )

        # Loop over combinations of fixing lower/upper bounds of other variables
        for other_symb_dict in list_of_other_symb_dicts:

            # Retrieve coordinates of specific hypercube enclosure line
            symb_vals = get_enclosure_line_coords_ndim(
                tm,
                symb,
                other_symb_dict,
                resolution=resolution,
            )
            coords_list.append(symb_vals)

    return cp1, coords_list


def plot_n_to_2_solution_enclosure(
    tm1,
    tm2,
    ax=None,
    resolution=100,
    color="k",
    show_center_point=False,
    plot_remainder_bound=True,
    verbose=False,
):
    """
    Plot the solution enclosure between two Taylor models in 2D.

    Args:
        tm1 (taylor_model): The first Taylor model.
        tm2 (taylor_model): The second Taylor model.
        ax: The axes to plot on.
        resolution (int): The number of points to use for plotting.
        color (str): The color of the plot.
        show_center_point (bool): Whether to show the center point.
        plot_remainder_bound (bool): Whether to plot the remainder bound.
        verbose (bool): Whether to print verbose output.

    Returns:
        ax: The updated axes with the plotted solution enclosure.
    """

    # Expand symbol set so taylor models include same set
    symbs1 = set(tm1.tpol.symbol_set)
    symbs2 = set(tm2.tpol.symbol_set)
    not_in_1 = symbs2 - symbs1
    not_in_2 = symbs1 - symbs2

    new_domains = {var: tm2.domain[var] for var in not_in_1}
    new_exp_points = {var: tm2.exp_point[var] for var in not_in_1}
    tm1.extend_symbol_set([var for var in not_in_1], new_exp_points, new_domains)

    new_domains = {var: tm1.domain[var] for var in not_in_2}
    new_exp_points = {var: tm1.exp_point[var] for var in not_in_2}
    tm2.extend_symbol_set([var for var in not_in_2], new_exp_points, new_domains)

    # Get coordinate lists for all combinations of hypercubes
    cp1, coords_list1 = get_solution_enclosure_ndim(tm1, resolution=resolution)
    cp2, coords_list2 = get_solution_enclosure_ndim(tm2, resolution=resolution)

    assert len(coords_list2) == len(coords_list1)

    # Loop through all the points
    for ls1, ls2 in zip(coords_list1, coords_list2):

        # Plot solution enclosure
        ax.plot(ls1, ls2, color=color, linewidth=0.5)

        if plot_remainder_bound:
            dx1 = tm1.rem_bound.upper - tm1.rem_bound.lower
            dx2 = tm2.rem_bound.upper - tm2.rem_bound.lower

            # Return if rem bound is 0 or infty
            if dx1 == 0 and dx2 == 0:
                if verbose:
                    print("Remainder bound is 0")
                continue
            if np.isinf(dx1) or np.isinf(dx2):
                if verbose:
                    print("Remainder bound is infinity.")
                continue

            # Get points from all the corners of the squares
            points = []
            for it2, (x1, x2) in enumerate(zip(ls1, ls2)):

                if it2 == 0 or it2 == len(ls1) - 1:
                    lower_x1 = x1 + tm1.rem_bound.lower
                    lower_x2 = x2 + tm2.rem_bound.lower
                    rect = Rectangle(
                        (lower_x1, lower_x2),
                        dx1,
                        dx2,
                        edgecolor=color,
                        facecolor="none",
                        alpha=1,
                        linestyle="--",
                        linewidth=0.5,
                    )
                    ax.add_patch(rect)

                # Corners of the rectangle
                x1 = x1 - dx1 / 2
                x2 = x2 - dx2 / 2
                corners = [
                    (x1, x2),
                    (x1 + dx1, x2),
                    (x1 + dx1, x2 + dx2),
                    (x1, x2 + dx2),
                ]

                # Add the corners to the list of points
                points.append(corners)

            points = np.array(points)
            for i in range(4):
                for j in range(len(points) - 1):
                    ax.plot(
                        points[j : j + 2, i, 0],
                        points[j : j + 2, i, 1],
                        color,
                        lw=0.5,
                        linestyle="--",
                        alpha=1,
                    )

    if show_center_point:
        ax.scatter(cp1, cp2, color=color)

    return ax


def plot_n_to_3_solution_enclosure(
    tm1,
    tm2,
    tm3,
    ax=None,
    resolution=100,
    color="k",
    show_center_point=False,
    plot_remainder_bound=True,
    verbose=False,
):
    """
    Plot the solution enclosure between two Taylor models in 3D.

    Args:
        tm1 (taylor_model): The first Taylor model.
        tm2 (taylor_model): The second Taylor model.
        ax: The axes to plot on.
        resolution (int): The number of points to use for plotting.
        color (str): The color of the plot.
        show_center_point (bool): Whether to show the center point.
        plot_remainder_bound (bool): Whether to plot the remainder bound.
        verbose (bool): Whether to print verbose output.

    Returns:
        ax: The updated axes with the plotted solution enclosure.
    """

    # Expand symbol set so taylor models include same set
    symbs1 = set(tm1.tpol.symbol_set)
    symbs2 = set(tm2.tpol.symbol_set)
    symbs3 = set(tm3.tpol.symbol_set)
    not_in_1 = (symbs2 | symbs3) - symbs1
    not_in_2 = (symbs1 | symbs3) - symbs2
    not_in_3 = (symbs1 | symbs2) - symbs3

    tm_23_domain = tm2.domain | tm3.domain
    new_domains = {var: tm_23_domain[var] for var in not_in_1}
    tm_23_exp_point = tm2.exp_point | tm3.exp_point
    new_exp_points = {var: tm_23_exp_point[var] for var in not_in_1}
    tm1.extend_symbol_set([var for var in not_in_1], new_exp_points, new_domains)

    tm_13_domain = tm1.domain | tm3.domain
    new_domains = {var: tm_13_domain[var] for var in not_in_2}
    tm_13_exp_point = tm1.exp_point | tm3.exp_point
    new_exp_points = {var: tm_13_exp_point[var] for var in not_in_2}
    tm2.extend_symbol_set([var for var in not_in_2], new_exp_points, new_domains)

    tm_21_domain = tm2.domain | tm1.domain
    new_domains = {var: tm_21_domain[var] for var in not_in_3}
    tm_21_exp_point = tm2.exp_point | tm1.exp_point
    new_exp_points = {var: tm_21_exp_point[var] for var in not_in_3}
    tm3.extend_symbol_set([var for var in not_in_3], new_exp_points, new_domains)

    # Get coordinate lists for all combinations of hypercubes
    cp1, coords_list1 = get_solution_enclosure_ndim(tm1, resolution=resolution)
    cp2, coords_list2 = get_solution_enclosure_ndim(tm2, resolution=resolution)
    cp3, coords_list3 = get_solution_enclosure_ndim(tm3, resolution=resolution)

    assert len(coords_list2) == len(coords_list1)
    assert len(coords_list3) == len(coords_list1)

    # Loop through all the points
    for ls1, ls2, ls3 in zip(coords_list1, coords_list2, coords_list3):

        # Plot solution enclosure
        ax.plot(ls1, ls2, ls3, color=color, linewidth=0.5)

        if plot_remainder_bound:
            dx1 = tm1.rem_bound.upper - tm1.rem_bound.lower
            dx2 = tm2.rem_bound.upper - tm2.rem_bound.lower
            dx3 = tm3.rem_bound.upper - tm3.rem_bound.lower

            # Return if rem bound is 0 or infty
            if dx1 == 0 and dx2 == 0 and dx3 == 0:
                if verbose:
                    print("Remainder bound is 0")
                continue
            if np.isinf(dx1) or np.isinf(dx2) or np.isinf(dx3):
                if verbose:
                    print("Remainder bound is infinity.")
                continue

            # Get points from all the corners of the squares
            points = []
            for it2, (x1, x2, x3) in enumerate(zip(ls1, ls2, ls3)):
                # Corners of the rectangle
                lower_x1 = x1 + tm1.rem_bound.lower
                lower_x2 = x2 + tm2.rem_bound.lower
                lower_x3 = x3 + tm3.rem_bound.lower
                vertices = [
                    (lower_x1, lower_x2, lower_x3),
                    (lower_x1 + dx1, lower_x2, lower_x3),
                    (lower_x1 + dx1, lower_x2 + dx2, lower_x3),
                    (lower_x1, lower_x2 + dx2, lower_x3),
                    (lower_x1, lower_x2, lower_x3 + dx3),
                    (lower_x1 + dx1, lower_x2, lower_x3 + dx3),
                    (lower_x1 + dx1, lower_x2 + dx2, lower_x3 + dx3),
                    (lower_x1, lower_x2 + dx2, lower_x3 + dx3),
                ]

                # Add the corners to the list of points
                points.append(vertices)

            ### Face collection ###
            faces = [
                (0,1,2,3),  # bottom
                (4,5,6,7),  # top
                (0,1,5,4),  # side
                (1,2,6,5),
                (2,3,7,6),
                (3,0,4,7)
            ]
            faces_list = []
            for it, cube in enumerate(points):
                #plot the cubes
                for face in faces:
                    faces_list.append([cube[face[0]], cube[face[1]], cube[face[2]], cube[face[3]]])
            
            # Optional: connect corresponding vertices as thin faces (quads) between cubes
                if it < len(points) - 1:
                    cube_a = points[it]
                    cube_b = points[it+1]
                    for v0, v1 in [(0,1),(1,2),(2,3),(3,0), (4,5),(5,6),(6,7),(7,4), (0,4),(1,5),(2,6),(3,7)]:
                        # create a quad connecting the edge from cube_a to cube_b
                        faces_list.append([cube_a[v0], cube_a[v1], cube_b[v1], cube_b[v0]])

            
            # Create Poly3DCollection
            pc = Poly3DCollection(faces_list, facecolors='cyan', edgecolors='k', linewidths=0.0005,
                                  alpha=0.005)
            ax.add_collection3d(pc)

    if show_center_point:
        ax.scatter(cp1, cp2, cp3, color=color)

    return ax


def sample_on_square_boundary(p, domain_size, n_samples_mc, margin=0.0, rng=None):
    """
    Samples points on the square boundary of a given domain. This function is used during Monte Carlo simulations to sample points more effectively.

    Args:
        p (tuple): The center point of the square (px, py).
        domain_size (float): The size of the domain.
        n_samples_mc (int): The number of samples to generate.
        margin (float): The margin from the edge of the domain.
        rng (np.random.Generator): Optional random number generator.

    Returns:
        np.ndarray: An array of shape (2, n_samples_mc) containing the sampled points.
    """
    # Specify random number generator
    if rng is None:
        rng = np.random.default_rng()

    px, py = p
    half = domain_size / 2 - margin  # Pull in slightly from the actual edge

    x_min, x_max = px - half, px + half
    y_min, y_max = py - half, py + half

    n_side = n_samples_mc // 4

    samples = []
    # Bottom side: y = y_min
    xs = rng.uniform(x_min, x_max, n_side)
    samples.extend(zip(xs, [y_min] * n_side))

    # Top side: y = y_max
    xs = rng.uniform(x_min, x_max, n_side)
    samples.extend(zip(xs, [y_max] * n_side))

    # Left side: x = x_min
    ys = rng.uniform(y_min, y_max, n_side)
    samples.extend(zip([x_min] * n_side, ys))

    # Right side: x = x_max
    ys = rng.uniform(y_min, y_max, n_side)
    samples.extend(zip([x_max] * n_side, ys))

    return np.array(samples).T  # Shape: (2, n_samples_mc)
