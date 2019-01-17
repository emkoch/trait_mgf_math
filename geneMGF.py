from copy import deepcopy
import itertools
import sympy
import re

def constGeneMGF(lineages):
    n_lineages = len(lineages)
    if n_lineages == 1:
        return sympy.Integer(1)
    dummy_vars = sympy.symbols(['s_' + '.'.join([str(tip) for tip in sorted(lineage)])
                                for lineage in lineages])
    mgf_factor = sympy.binomial(n_lineages, 2)
    for ss in dummy_vars:
        mgf_factor -= ss
    coal_events = [cc for cc in itertools.combinations(lineages, 2)]
    result = sympy.Integer(0)
    for coal_event in coal_events:
        new_lineages = lineages[:]
        for lineage in coal_event:
            new_lineages.remove(lineage)
        new_lineages.append(sorted(coal_event[0] + coal_event[1]))
        result += mgf_factor**-1*constGeneMGF(new_lineages)
    return result

def partition(collection):
    if len(collection) == 1:
        yield collection
        return
    first = collection[0]
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartitions subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [first  + subset]  + smaller[n+1:]
        # put `first` in its own subset
        yield [first] + smaller

def deme_partition(collection, n_demes):
    if len(collection) == 1:
        for deme_index in range(n_demes):
            result = []
            for ii in range(deme_index):
                result.append([])
            result.append(collection)
            for ii in range(deme_index + 1, n_demes):
                result.append([])
            yield result
        return
    first = collection[0]
    for smaller_deme_part in deme_partition(collection[1:], n_demes):
        for deme in range(n_demes):
            result = deepcopy(smaller_deme_part)
            result[deme] = sorted(result[deme] + [first])
            yield result

def deme_part_to_symbol(deme_part):
    def lin_to_str(lineage):
        return '(' + '.'.join(sorted([str(indiv) for indiv in lineage])) + ')'
    result = 'phi_'
    result += '.'.join(['[' + ';'.join(sorted([lin_to_str(lineage) for lineage in deme])) + ']'
                        for deme in deme_part])
    return sympy.symbols(result)

def deme_symbol_to_part(deme_symbol):
    parts = deme_symbol.split("phi_")[1].split("[")[1:]
    p = re.compile("\((.*?)\)")
    return [[lin.split(".") for lin in p.findall(part)] for part in parts]

def deme_part_symbol_after_coal(deme_part, coal_pair):
    result = deepcopy(deme_part)
    for deme_idx, deme in enumerate(result):
        if coal_pair[0] in deme and coal_pair[1] in deme:
            result[deme_idx].remove(coal_pair[0])
            result[deme_idx].remove(coal_pair[1])
            result[deme_idx].append(sorted(coal_pair[0] + coal_pair[1]))
    return deme_part_to_symbol(result)

def deme_part_symbol_after_migr(deme_part, lineage, deme_from, deme_to):
    result = deepcopy(deme_part)
    result[deme_from].remove(lineage)
    result[deme_to].append(lineage)
    return deme_part_to_symbol(result)

class geneMGF(object):
    def __init__(self, lineages, **kwargs):
        self.lineages = lineages
        self.num_lineages = len(lineages)
        self.mgf = self.create_mgf_expr(**kwargs)
        self.mgf_type = 'standard'

    def create_mgf_expr(self):
        """ Generates the mgf, only should be called by __init__ """
        if self.num_lineages == 1:
            return sympy.Integer(1)
        dummy_vars = sympy.symbols(['s_' + '.'.join([str(tip) for tip in sorted(lineage)])
                                    for lineage in self.lineages])
        mgf_factor = sympy.binomial(self.num_lineages, 2)
        for ss in dummy_vars:
            mgf_factor -= ss
        coal_events = [cc for cc in itertools.combinations(self.lineages, 2)]
        result = sympy.Integer(0)
        for coal_event in coal_events:
            new_lineages = self.lineages[:]
            for lineage in coal_event:
                new_lineages.remove(lineage)
            new_lineages.append(sorted(coal_event[0] + coal_event[1]))
            result += mgf_factor**-1*constGeneMGF(new_lineages)
        return result

    def calc_gene_mom(self, branches, orders):
        """ Calculate a moment of the genealogy
        Arguments:
        branches -- which branches should be included in the moment calculation
                            ex: [1]     ex: ['1.2', 3]
        orders   -- the orders corresponding to the branches in the moment calculation
                    list of integers
        """
        if not isinstance(branches, list):
            branches = list(branches)
        if not isinstance(orders, list):
            orders = list(orders)
        dummy_vars = sympy.symbols(['s_' + str(branch) for branch in branches])
        d_mgf = self.mgf
        for branch_idx in range(len(branches)):
            d_mgf = sympy.diff(d_mgf, dummy_vars[branch_idx], orders[branch_idx])
        all_dummy_vars = sympy.symbols([str(fs) for fs in d_mgf.free_symbols if 's_' in str(fs)])
        return d_mgf.subs([(dummy_var, 0) for dummy_var in all_dummy_vars])

    def gene_mgf_to_trait(self, moment_order):
        """ subs into the gene mgf to get the trait mgf
        Arguments:
        moment_order -- determines the order to which the mutational mgf is expanded
        """
        branches = [indv_combo for branch_size in range(1, self.num_lineages)
                    for indv_combo in itertools.combinations(self.lineages, branch_size)]
        result = self.mgf
        for branch in branches:
            branch_term = sympy.symbols('s_' + '.'.join([str(tip[0]) for tip in branch]))
            branch_sum = sympy.S(0)
            for indiv in branch:
                branch_sum += sympy.symbols('k_' + str(indiv[0]))
            mut_term = sympy.S(0)
            for mut_order in range(1, moment_order + 1):
                mut_term += (sympy.symbols('m_' + str(mut_order))/sympy.factorial(mut_order)*
                             branch_sum**mut_order)
            result = result.subs(branch_term, sympy.symbols('theta')*mut_term)
        return result

class geneMGF_struct_slow(geneMGF):
    def __init__(self, lineages, **kwargs):
        geneMGF.__init__(self, lineages, **kwargs)
        self.mgf_type = 'structured'

    def create_mgf_expr(self, coal_rates, migration_rates, initial_config):
        """ Generates the mgf, only should be called by __init__
        Arguments:
        coal_rates      -- a list giving the coalescent rate in each deme
                           ex: [1, 1]
        migration_rates -- a matrix giving the migration rates from each deme to each other deme
                           ex: sympy.Matrix([[0, 1], [1, 0]])
        initial_config  -- a list giving the starting configuration of the lineages
                           ex: [ [[1], [2]], [[3]] ]
        """
        Eqns = []
        num_demes = len(coal_rates)
        # Get all partitions of starting lineages aka all coalescent states
        indv_partitions = list(partition(self.lineages))
        indv_partitions.remove([[individual for lineage in self.lineages
                                 for individual in lineage]])
        deme_partitions = [deme_part for indv_partition in indv_partitions
                           for deme_part in deme_partition(indv_partition, num_demes)]
        print("There are " + str(len(deme_partitions)) + " possible states." )
        print("Setting up a system of " + str(len(deme_partitions)) + " equations...")
        deme_part_symbols = []
        all_symbols = []
        for deme_part in deme_partitions:
            print("Working with partition: " + str(deme_part))
            deme_part_symbol = deme_part_to_symbol(deme_part)
            deme_part_symbols.append(deme_part_symbol)
            all_symbols.append(deme_part_symbol)
            # Make the factor for this partition
            print("Creating the factor for this partition")
            deme_part_factor = sympy.Integer(0)
            for deme_idx, deme in enumerate(deme_part):
                n_deme_lineages = len(deme)
                deme_part_factor += sympy.binomial(n_deme_lineages, 2)*coal_rates[deme_idx]
                for other_deme_idx in range(num_demes):
                    if deme_idx != other_deme_idx:
                        deme_part_factor += n_deme_lineages*migration_rates[deme_idx,
                                                                            other_deme_idx]
                for lineage in deme:
                    deme_part_factor -= sympy.symbols('s_' +
                                                      '.'.join([str(tip)
                                                                for tip in sorted(lineage)]))
            # Add the left side in (sign is negative, see eq 8 Lohse et al.)
            deme_part_eqn = -deme_part_factor*deme_part_symbol
            # sympy.pprint(deme_part_eqn)
            # Add the transfer rates to coalescent states
            print("Adding rates to other coal states")
            for deme_idx, deme in enumerate(deme_part):
                for pair in itertools.combinations(deme, 2):
                    # If only two lineages left, the next coal event reduces us to one
                    if sum([len(deme_lins) for deme_lins in deme_part]) == 2:
                        symbol_after_coal = sympy.Integer(1)
                    else:
                        symbol_after_coal = deme_part_symbol_after_coal(deme_part, pair)
                        all_symbols.append(symbol_after_coal)
                    deme_part_eqn += coal_rates[deme_idx]*symbol_after_coal
            # Add the transfer rates for migration events
            print("Adding rates to other migration states")
            for deme_idx, deme in enumerate(deme_part):
                # For each lineage in the deme add a term for the probability of that lineage
                # immigrating to each other deme
                for lin in deme:
                    for deme_to_idx in range(num_demes):
                        deme_part_eqn += (migration_rates[deme_idx, deme_to_idx]*
                                          deme_part_symbol_after_migr(deme_part, lin,
                                                                      deme_idx, deme_to_idx))
                        all_symbols.append(deme_part_symbol_after_migr(deme_part, lin,
                                                                      deme_idx, deme_to_idx))
            print("Done! Adding this to the list of equations")
            Eqns.append(deme_part_eqn)
        print("The terms being solved for are:")
        for deme_part_symbol in deme_part_symbols:
            sympy.pprint(deme_part_symbol)
        all_symbols = set(all_symbols)
        print("There are " + str(len(all_symbols)) + " terms in the system of equation")
        for symbol in all_symbols:
            sympy.pprint(symbol)
        starting_part = deme_part_to_symbol(initial_config)
        starting_idx = deme_part_symbols.index(starting_part)
        print('Solving system of equations...')
        mgf_solution = sympy.linsolve(Eqns, *deme_part_symbols)
        return list(mgf_solution)[0][starting_idx]