# this file was created by Piotr Gawron it contains large portions of code from
# https://github.com/lanl-ansi/bqpsolvers/ repository

import warnings
from gurobipy import Model, GRB
import dimod


class GurobiSampler(dimod.Sampler):
    properties = None
    parameters = None
    runnable = None

    def __init__(self):
        self.parameters = {
            'method': []
        }

        self.properties = {}

    def sample(self, bqm:dimod.BinaryQuadraticModel, method="miqp", num_reads=1, gurobi_params_kw=None):
        assert(method in ["mip", "miqp"])
        bqm_bin = bqm.change_vartype(vartype=dimod.BINARY, inplace=False)
        variable_ids = frozenset(bqm_bin.variables)
        variable_product_ids = frozenset(bqm_bin.quadratic)

        m = Model()

        gurobi_params = {
            "OutputFlag": 0,
            "TimeLimit": 60,
            "Threads": 12,
            "Cuts": 1,
            "MIPFocus": 2,
            "PoolSearchMode": 2,
            "PoolSolutions": num_reads
        }

        if gurobi_params_kw is None:
            gurobi_params_kw = {}

        gurobi_params.update(gurobi_params_kw)

        for param, value in gurobi_params.items():
            m.setParam(param, value)

        variable_lookup = {}
        for vid in variable_ids:
            variable_lookup[vid] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY, name="var_{}".format(vid))
        if method == "mip":
            for pair in variable_product_ids:
                variable_lookup[pair] = m.addVar(lb=0, ub=1, vtype=GRB.BINARY,
                                                 name="link_{}_{}".format(str(pair[0]), str(pair[1])))
        m.update()
        if method == "mip":
            for i, j in variable_product_ids:
                m.addConstr(variable_lookup[(i, j)] >= variable_lookup[i] + variable_lookup[j] - 1)
                m.addConstr(variable_lookup[(i, j)] <= variable_lookup[i])
                m.addConstr(variable_lookup[(i, j)] <= variable_lookup[j])

        bqm_ising = bqm.change_vartype(vartype=dimod.SPIN, inplace=False)
        if len(bqm_ising.linear) <= 0 or all(bqm_ising.linear[lt] == 0.0 for lt in bqm_ising.linear):
            warnings.warn('detected spin symmetry, adding symmetry breaking constraint')
            v1 = variable_ids[0]
            m.addConstr(variable_lookup[v1] == 0)

        obj = 0.0
        for lt in bqm_bin.linear:
            obj += bqm_bin.linear[lt] * variable_lookup[lt]

        if method == "mip":
            for qt in bqm_bin.quadratic:
                obj += bqm_bin.quadratic[qt] * variable_lookup[qt]
        elif method == "miqp":
            for qt in bqm_bin.quadratic:
                i = qt[0]
                j = qt[1]
                obj += bqm_bin.quadratic[qt] * variable_lookup[i] * variable_lookup[j]

        m.setObjective(obj, GRB.MINIMIZE)
        m.update()

        m.optimize()

        energies = []
        samples = []
        for i in range(m.SolCount):
            m.Params.SolutionNumber = i  # set solution numbers
            energy = m.ObjVal + bqm_bin.offset
            sample = {k: int(variable_lookup[k].X) for k in variable_ids}
            energies.append(energy)
            samples.append(sample)

        ss = dimod.SampleSet.from_samples(samples, vartype=bqm.BINARY, energy=energies, aggregate_samples=True)
        return ss.change_vartype(bqm.vartype)


if __name__ == "__main__":
    import random
    import itertools
    num_vars = 10
    linear = {k: v for k, v in enumerate(random.random() for _ in range(num_vars))}
    quadratic = {(k1, k2): random.random() for k1, k2 in itertools.product(range(num_vars), range(num_vars)) if k1 != k2}

    bqm = dimod.BinaryQuadraticModel(linear,
                                     quadratic,
                                     offset=3.0,
                                     vartype=dimod.SPIN)

    sampler = GurobiSampler()
    print(sampler.sample(bqm, method="mip", num_reads=20, gurobi_params_kw={"TimeLimit": 1}))
    print(sampler.sample(bqm, method="miqp", num_reads=20, gurobi_params_kw={"TimeLimit": 1}))

    if num_vars <= 10:
        sampler = dimod.ExactSolver()
        print(sampler.sample(bqm))
