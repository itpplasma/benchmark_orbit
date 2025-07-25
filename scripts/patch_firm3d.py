#!/usr/bin/env python3
"""
Patch firm3d to work with newer booz_xform that uses toroidal_flux instead of chi.
"""

import os

# Path to the boozermagneticfield.py file
boozer_file = "codes/firm3d/src/simsopt/field/boozermagneticfield.py"

# Read the file
with open(boozer_file, 'r') as f:
    content = f.read()

# Replace chi with toroidal_flux
content = content.replace("self.bx.chi", "self.bx.toroidal_flux")

# Write back
with open(boozer_file, 'w') as f:
    f.write(content)

print(f"Patched {boozer_file} to use toroidal_flux instead of chi")