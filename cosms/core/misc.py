#!/usr/bin/env python3
# coding: utf-8
"""misc
"""
import math

MAD_FACTOR = 1.4826


def median(values, sorting = True):
    """Calculate median of list of values
    
    :param list values: List of floating point numbers
    :param bool sorting: Values will be sorted if True
    
    :return median_value: Median
    :rtype: float
    """
    if sorting:
        values.sort()
    n = len(values)
    if n%2 == 1:
        median_value = values[math.floor(n/2)]
    else:
        median_value = sum(values[int((n/2)-1):int((n/2)+1)]) / 2.0
    return median_value


def median_absolute_deviation(values, sorting = True):
    """Median absolute deviation (MAD)
    
    :param list values: List of floating point numbers
    :param bool sorting: Values will be sorted if True
    
    :return mad: Median absolute deviation
    :rtype: float
    """
    if sorting:
        values.sort()
    median_value = median(values, sorting = False)
    mad = median(sorted([abs(i - median_value) for i in values]), 
                 sorting = False)
    return mad


def min_max_median_mean(values, sorting = True):
    """Basic descriptive statistics

    :param list values: List of values
    :param bool sorting: Values will be sorted if True

    :return min_value:
    :return max_value:
    :return median_value:
    :return arithmetic_mean:
    :return standard_deviation:
    """
    if sorting:
        values.sort()
    min_value = values[0]
    max_value = values[-1]
    n = len(values)
    total = sum(values)
    arithmetic_mean = total / n
    median_value = median(values, sorting = False)
    standard_deviation = math.sqrt(sum([(i - arithmetic_mean)**2 
                                        for i in values]) / n)
    return min_value, max_value, median_value, arithmetic_mean, standard_deviation


def std_dev_from_mad(mad):
    """Estimate standard deviation from median absolute deviation (MAD)
    
    :param float mad: Median absolute deviation
    
    :return std_dev: Standard deviation
    :rtype: int
    """
    std_dev = mad * MAD_FACTOR
    return std_dev


class LazyFunction(object):
    """This replaces a function by its returned value so that the original 
    function is calculated only once. It also prevents asking whether or not
    the function has already been called before.
    """
    def __init__(self, func):
        self.func = func
        

    def __get__(self, obj, typ):
        if obj is None:
            return self
        value = self.func(obj)
        setattr(obj, self.func.__name__, value)
        return value
    