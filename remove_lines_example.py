#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 13:25:03 2022

@author: boeglinw
"""

import matplotlib.pyplot as plt




line1 = [2, 4, 8]
line2 = [3, 6, 12]

fig, ax = plt.subplots(1,1)

ax.plot(line1, label = 'line1')
ax.plot(line2, linestyle='dashed', label = 'line2')
#(line_2.pop(0)).remove()



plt.legend()

ax.get_legend().remove()

# line_2.pop().remove()

plt.margins(0.2)
plt.title("With extra lines")
plt.show()

"""
plt.plot(line1)
l = line_2.pop(0)
l.remove()
plt.margins(0.2)
plt.title("Removed extra lines")
"""
plt.show()