import gqcpy
import prototypes.vector


def magnetic_field_along_axis(B, n):
    """Create a HomogeneousMagneticField with strength B along the axis n."""
    return gqcpy.HomogeneousMagneticField(B * n)


def magnetic_field_along_random_axis(B):
    """Create a HomogeneousMagneticField with strength B along a random axis."""
    n = prototypes.vector.create_random_unit_vector()
    return magnetic_field_along_axis(B, n)


def create_magnetic_fields(B_values):
    magnetic_fields = []
    for B in B_values:
        magnetic_fields.append(gqcpy.HomogeneousMagneticField([0,0,-B]))

    return magnetic_fields
