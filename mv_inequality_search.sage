def V(A: Polyhedron, B: Polyhedron):
    """
    Compute the mixed volume of two polyhedra A and B in R^2.
    """
    return 1/2 * ((A+B).volume() - A.volume() - B.volume())

# Generate polynomial inequality terms
def generate_inequality_terms(n: int):
    """
    Generate all possible terms of polynomial inequality terms
    of the form V(a,b)V(c,d) where a neq b neq c \neq d are elements
    of an n-body set.
    """
    bodies = Set(range(n))
    term_list = []
    for mv_1 in Subsets(bodies, 2):
        if n - mv_1[0] < 2:
            # Not enough space for another mv factor
            continue
        
        # Get subset of remaining bodies so that terms aren't repeated
        remaining = Set(range(mv_1[0]+1, n)) - [mv_1[1]]

        for mv_2 in Subsets(remaining, 2):
            term_list.append(Set([mv_1, mv_2]))
    
    assert len(term_list) == binomial(n, 2) * binomial(n-2, 2) / 2

    return Set(term_list)

def compute_inequality_terms(terms: Set, bodies: Set):
    """
    Compute the values of polynomial inequality terms of the form V(a,b)V(c,d)
    where a neq b neq c neq d provided a set of bodies.
    """
    return tuple([
        V(bodies[term[0][0]], bodies[term[0][1]]) * V(bodies[term[1][0]], bodies[term[1][1]])
        for term in terms
    ])

body_types = Set([
    Polyhedron([(0,0), (1,0)]),
    Polyhedron([(0,0), (0,1)]) 
])

def compute_rays(n: int):
    """
    Find the rays in term value space produced by n bodies when the terms follow
    the form V(a,b)V(c,d) where a neq b neq c neq d.
    """
    terms = generate_inequality_terms(n)
    body_space = cartesian_product([body_types] * n)
    return list(filter(
        lambda x: sum(x) != 0, 
        map(
            lambda s: compute_inequality_terms(terms, s),
            body_space
        )
    ))

def compute_convex_hull(n: int):
    """
    Find convex hull of rays in term value space produced by n bodies when the terms follow
    the form V(a,b)V(c,d) where a neq b neq c neq d.
    """
    rays = compute_rays(n)
    pmv = Polyhedron(rays=rays)
    return pmv