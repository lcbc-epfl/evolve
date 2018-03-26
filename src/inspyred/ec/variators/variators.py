"""
    ================
    :mod:`variators`
    ================
    
    .. Copyright 2012 Aaron Garrett

    .. This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

    .. This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

    .. You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>.
       
    .. module:: variators
    .. moduleauthor:: Aaron Garrett <aaron.lee.garrett@gmail.com>
"""


def default_variation(random, candidates, args):
    """Return the set of candidates without variation.

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments
    
    """
    return candidates

