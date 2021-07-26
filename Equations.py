import math

class Spring_Eq:

    def __init__(self, k0: float, dk: float, w: float, theta: float, mass=1):
        self.k0 = k0
        self.dk = dk
        #self.w = math.sqrt(self.k0 / mass)
        self.w = w
        self.theta = theta

    def change_theta(self, theta: float):
        """returns a new spring equation with the new theta value"""
        return Spring_Eq(self.k0, self.dk, theta)

    def calc(self, time):
        """calculates the spring constant value at a specified time from the spring equation, and returns
        that float number"""
        return (self.k0 + (self.dk * math.sin(((self.w * time) + self.theta))))

    def __str__(self):
        return(f"k = {self.k0} + ({self.dk} * sin(({self.w} * t) + {self.theta}))")


def get_spring_equations(k0: float, dk: float, w: float, thetas: list):
    """returns a list of spring equations with the different thetas suplied in the theta list"""
    equations = []
    for theta in thetas:
        eq = Spring_Eq(k0, dk, w, theta)
        equations.append(eq)
    return equations

class Spring_Eq_2:
    """
    Spring equations based off of the formula
        k = k0 * [1 + dk*cos((2pi / N)i - w*t)
    """

    def __init__(self, k0: float, dk: float, w: float, i: int, N: int, mass=1):
        self.k0 = k0
        self.dk = dk
        #self.w = math.sqrt(self.k0 / mass)
        self.w = w
        self.i = i
        self.N = N


    def calc(self, time):
        """calculates the spring constant value at a specified time from the spring equation, and returns
        that float number"""
        return (self.k0 * (1 + (self.dk * math.cos((2 * math.pi * self.i / self.N) - (self.w * time)))))

    def __str__(self):
        return(f"k = {self.k0}[1 + {self.dk} * cos((2pi / N)i - {self.w}*t)")


def get_spring_equations_2(k0: float, dk: float, w: float, N: int):
    """returns a list of spring equations"""
    equations = []
    for i in range(N):
        eq = Spring_Eq_2(k0, dk, w, i, N)
        equations.append(eq)
    return equations


class Force_Eq:

    def __init__(self, f0: float, df: float, w: float, theta: float):
        self.f0 = f0
        self.df = df
        self.w = w
        self.theta = theta

    def get_equations(self, f0: float, df: float, thetas: list):
        """returns a list of force equations with the different thetas suplied in the theta list"""
        equations = []
        for theta in thetas:
            eq = Spring_Eq(f0, df, theta)
            equations.append(eq)
        return equations

    def calc_mag(self, time):
        """calculates the force magnitude value at a specified time from the force equation, and returns
        that float number"""

        return (self.f0 + (self.df * math.sin((self.w * time) + self.theta)))

    def __str__(self):
        return(f"k = {self.f0} + ({self.df} * sin(({self.w} * t) + {self.theta}))")
