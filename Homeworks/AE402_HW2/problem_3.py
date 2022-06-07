"""
    problem_3_refactored.py created by Seb at 16:09 18/10/2021
    Path: HW2

    From the given information:

    (a) Compute and plot the position (non-dimensional units) of the five Lagrange points for the Sun-Jupiter
        system. Your plot must also include the positions of the primary and secondary bodies.

    (b) For the Jupiter-Europa system, compute the distance of L1 and L2 in km, with respect to Europa.

    (c) For the Sun-Apophis system, compute the distance of L1 in km, with respect to Apophis.

"""
from math import sqrt
import matplotlib.pyplot as plt, pandas as pd

from sympy import symbols, solve


def weighted_mu_value(mu_big: float, mu_small: float) -> float:
    """
    Simply return weighted values of both mu_s normalized to their sum
    :param mu_big:      mu of larger body
    :param mu_small:    mu of second body
    :return:            float weighted value
    """
    return mu_small / (mu_big + mu_small)


def absolute_values(lagrange_points: dict, one_du: float) -> dict:
    for point in lagrange_points:
        lagrange_points[point][0] = lagrange_points[point][0] * one_du
        lagrange_points[point][1] = lagrange_points[point][1] * one_du

    return lagrange_points


def calc_distance_between(point_1: list, point_2: list) -> float:
    # print("points: ", point_1, point_2)
    delta_x = abs(point_1[0] - point_2[0])
    delta_y = abs(point_1[1] - point_2[1])

    # print(delta_x, delta_y)
    distance = sqrt(delta_x ** 2 + delta_y ** 2)

    return distance


def lagrange_points(weighted_mu: float, filename: str) -> tuple:
    """
        Return lagrange points in distance units for the specified two-body system
        :param filename:            name of csv file to put points in
        :param weight_mu:           weight value of mu
        :return:                    list of lagrange points x-y pairs
    """
    mu = weighted_mu
    points = {}

    x1, x2, x3 = symbols('x1', real=True), symbols('x2', real=True), symbols('x3', real=True)

    l1_func = x1 - ((1 - mu) / (mu + x1) ** 2) + (mu / (x1 - 1 + mu) ** 2)
    l2_func = x2 - ((1 - mu) / (mu + x2) ** 2) - (mu / (x2 - 1 + mu) ** 2)
    l3_func = x3 + ((1 - mu) / (mu + x3) ** 2) + (mu / (x3 - 1 + mu) ** 2)

    # Applying the numerical solving solution
    l1_solution, l2_solution, l3_solution = solve(l1_func, x1), solve(l2_func, x2), solve(l3_func, x3)

    # Putting results into a list, all on the y=0 plane
    l1_coords, l2_coords, l3_coords = [l1_solution[0], 0], [l2_solution[0], 0], [l3_solution[0], 0]

    # appending all results to list, plus L4 and L5 which can be easily calculated
    points['L1'] = l1_coords
    points['L2'] = l2_coords
    points['L3'] = l3_coords
    points['L4'] = [0.5 - mu, sqrt(3) / 2]
    points['L5'] = [0.5 - mu, -sqrt(3) / 2]

    lagrange_dataframe = pd.DataFrame(columns=['Position', 'x', 'y'])
    for i in range(5):
        lagrange_dataframe.loc[i] = [f"L{i + 1}", points[f"L{i + 1}"][0], points[f"L{i + 1}"][1]]

    lagrange_dataframe.to_csv(filename)

    return points


def main():
    data = {"Mass": {
                "Sun": 1.989e30,
                "Earth": 5.972e24,
                "Moon": 7.347e22,
                "Jupiter": 1.898e27,
                "Europa": 4.799e22,
                "Apophis": 2.699e10
                },
            "Distances": {
                "Sun_Jupiter": 750e6,
                "Sun_Apophis": 149e6,
                "Earth_Moon": 384e3,
                "Jupiter_Europa": 670e3,
                }
            }

    # Lagrange points for the Sun-Jupiter System
    sun_jupiter_weighted_mu = weighted_mu_value(data['Mass']['Sun'], data['Mass']['Jupiter'])
    jupiter_europa_weighted_mu = weighted_mu_value(data['Mass']['Jupiter'], data['Mass']['Europa'])
    sun_apophis_weighted_mu = weighted_mu_value(data['Mass']['Sun'], data['Mass']['Apophis'])
    print(sun_apophis_weighted_mu)
    # sun_jupiter_lagrange_points = lagrange_points(sun_jupiter_weighted_mu, "problem_3a_lagrange.csv")
    # jupiter_europa_lagrange_points = lagrange_points(jupiter_europa_weighted_mu, "problem_3b_lagrange.csv")
    # sun_apophis_lagrange_points = lagrange_points(sun_apophis_weighted_mu, "problem_3c_lagrange.csv")

    sun_jupiter_lagrange_df = pd.read_csv("problem_3a_lagrange.csv")
    jupiter_europa_lagrange_df = pd.read_csv("problem_3b_lagrange.csv")
    sun_apophis_lagrange_df = pd.read_csv("problem_3c_lagrange.csv")

    sun_jupiter_lagrange_points = {}
    for i, row in sun_jupiter_lagrange_df.iterrows():
        sun_jupiter_lagrange_points[row['Position']] = [row['x'], row['y']]

    jupiter_europa_lagrange_points = {}
    for i, row in jupiter_europa_lagrange_df.iterrows():
        jupiter_europa_lagrange_points[row['Position']] = [row['x'], row['y']]

    sun_apophis_lagrange_points = {}
    for i, row in sun_apophis_lagrange_df.iterrows():
        sun_apophis_lagrange_points[row['Position']] = [row['x'], row['y']]

    # For part b:

    # First convert all lagrange point dictionaries into absolute values
    jupiter_europa_absolute = absolute_values(jupiter_europa_lagrange_points,
                                              data['Distances']['Jupiter_Europa'])
    sun_apophis_absolute = absolute_values(sun_apophis_lagrange_points,
                                              data['Distances']['Sun_Apophis'])
    # Calculate distances for Jupiter-Europa
    europa_to_l1 = calc_distance_between(jupiter_europa_absolute['L1'],
                                         [data['Distances']['Jupiter_Europa'], 0])
    europa_to_l2 = calc_distance_between(jupiter_europa_absolute['L2'],
                                         [data['Distances']['Jupiter_Europa'], 0])

    # Calculate distances for Sun-Aphosis
    apophis_to_l1 = calc_distance_between(sun_apophis_absolute['L1'],
                                         [data['Distances']['Sun_Apophis'], 0])

    intro_string = f"HW2 Problem 3 - Lagrange points in different Systems.\n" \
                   f"------------------------------------------------------\n"

    # Part 3a
    pt_a_string = f"Part 3a: Sun-Jupiter Lagrange points:\n"
    for point in sun_jupiter_lagrange_points:
        pt_a_string += f"{point}: [{sun_jupiter_lagrange_points[point][0]:.3f}, {sun_jupiter_lagrange_points[point][1]:.3f}]\n"
    pt_a_string += f"------------------------------------------------------\n"

    # Part 3b string
    pt_b_string = f"Part 3b: Jupiter-Europa Lagrange points:\n"
    for point in jupiter_europa_lagrange_points:
        pt_b_string += f"{point}: [{jupiter_europa_lagrange_points[point][0]:.3f}, {jupiter_europa_lagrange_points[point][1]:.3f}]\n"
    pt_b_string += f"------------------------------------------------------\n" \
                   f"Mu value of {jupiter_europa_weighted_mu:.10f}\n" \
                   f"Distance from Europa to Jupiter-Europa L1: {europa_to_l1:.3f} km\n" \
                   f"Distance from Europa to Jupiter-Europa L2: {europa_to_l2:.3f} km\n" \
                   f"------------------------------------------------------\n"

    # Part 3c string
    pt_c_string = f"Part 3c: Sun-Apophis Lagrange points:\n"
    for point in sun_apophis_lagrange_points:
        pt_c_string += f"{point}: [{sun_apophis_lagrange_points[point][0]:.3f}, {sun_apophis_lagrange_points[point][1]:.3f}]\n"
    pt_c_string += f"------------------------------------------------------\n" \
                   f"Distance from Apophis to Sun-Apophis L1: {apophis_to_l1:,.3f} km\n"

    # Writing output to file for future reference
    with open('output/problem_3_output.txt', 'w') as output:
        output.write(intro_string)
        output.write(pt_a_string)
        output.write(pt_b_string)
        output.write(pt_c_string)

    print(intro_string + pt_a_string + pt_b_string + pt_c_string)

    # Setting up plot
    fig1, ax1 = plt.subplots()

    # Plotting Lagrange points
    for i, row in sun_jupiter_lagrange_df.iterrows():
        ax1.scatter(row['x'], row['y'], c='k', s=5)

    # Plotting Sun and Jupiter Positions
    ax1.scatter([0], [0], s=600, c='yellow', label="Sun")
    ax1.scatter([1], [0], s=100, c='orange', label="Jupiter")

    # Turning gridlines and legend on, setting ticks fontsize
    ax1.minorticks_on()
    ax1.grid(b=True, which='minor', color='k', linestyle='-', alpha=0.05)
    ax1.grid(b=True, which='major', color='k', linestyle='-', alpha=0.3)
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)

    # Making axis scales correctly
    leg = ax1.legend(scatterpoints=1, fontsize=10)
    leg.legendHandles[0]._sizes = [100]
    leg.legendHandles[1]._sizes = [100]

    # Setting x/y axis limits
    # ax1.set_xlim(-1.4, 1.4)
    # ax1.set_ylim(-0.6, 0.6)

    # Labels
    ax1.set_ylabel("Y-Axis DU units")
    ax1.set_xlabel("X-Axis DU units")
    ax1.set(title="Sun-Jupiter System Lagrange points\n(Planet size not to scale)")

    # saving to file
    fig1.savefig("problem_3a_lagrange_plot.png")
    # plt.show()

    return 0


if __name__ == '__main__':
    main()
