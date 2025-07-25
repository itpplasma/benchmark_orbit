#!/usr/bin/env python3
"""
Patch firm3d to work with newer booz_xform that doesn't have chi.
We'll calculate chi by integrating phip.
"""

import os

# Path to the boozermagneticfield.py file
boozer_file = "codes/firm3d/src/simsopt/field/boozermagneticfield.py"

# Read the file
with open(boozer_file, 'r') as f:
    content = f.read()

# Find the line that uses chi and replace with integrated phip
old_line = "        psip = self.bx.toroidal_flux / (2 * np.pi)"
new_lines = """        # Calculate chi by integrating phip
        s_half = 0.5 * (self.bx.s_b[1:] + self.bx.s_b[:-1])  # s values on half grid
        ds = np.diff(self.bx.s_b)
        chi = np.zeros(len(self.bx.s_b))
        chi[1:] = np.cumsum(self.bx.phip[1:] * ds)  # Integrate phip to get chi
        psip = chi / (2 * np.pi)"""

content = content.replace(old_line, new_lines)

# Write back
with open(boozer_file, 'w') as f:
    f.write(content)

print(f"Patched {boozer_file} to calculate chi from phip")