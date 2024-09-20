from graph import *
from graph_io import *
from graphviz import *
import random
from scipy.optimize import linprog, LinearConstraint, milp
import math
import numpy as np

def make_random_graph(n,m):
    """
    make random digraph with n nodes and m edges
    """
    if m<n:
        raise ValueError("Not enough edges!")
    G = MyGraph(True, n, True)
    Vlist = G.vertices
    for i in range(len(Vlist)):
        Vlist[i].label=i
    Vlist = G.vertices
    for v in Vlist:
        w=v
        while w==v:
            w = random.choice(Vlist)
        e=Edge(v,w)
        G.add_edge(e)
    l=m-n
    while l>0:
        v=random.choice(Vlist)
        w=random.choice(Vlist)
        if v != w:
            f=G.find_edge(v,w)
            f=list(f)
            if len(f) == 0 or (len(f) == 1 and f[0].tail==w):
                e=Edge(v,w)
                G.add_edge(e)
                l-=1
    return G

def make_cycle_list(G):
    """
    makes list of all cycles of a graph
    :param G: Graph
    :return: cycle_list
    """
    used_nodes = set()
    cycle_list = []
    for v in G.vertices:
        cycle_list += find_cycles(G,[v],used_nodes)
        used_nodes.add(v)
    return cycle_list

def generate_and_print_random_valid_sign_list(LL):
    #assign signs to cycles based on a valid pattern
    signlist = []
    n=len(G.edges)
    for e in G.edges:
        e.phiweight = random.choice([-1,1])*n**random.normalvariate(0,1)
    for L in LL:
        strr = ''
        tot = 0
        for i in range(len(L) - 1):
            strr += str(L[i].label) + '-'
            eset = G.find_edge(L[i], L[i+1])
            for e in eset:
                if e.tail == L[i]:
                    tot += e.phiweight
        eset = G.find_edge(L[-1], L[0])
        for e in eset:
            if e.tail == L[-1]:
                tot += e.phiweight
        strr += str(L[-1].label)
        if tot > 0:
            sgn = '+'
        elif tot <0:
            sgn = '-'
        else:
            raise ValueError('what are the odds that one of the cycles has weight 0?')
        signlist.append(sgn)
        print(strr + ' ', sgn)
    return signlist

def generate_and_print_random_sign_list(LL):
    signlist = []
    for L in LL:
        strr=''
        for i in range(len(L)-1):
            strr+=str(L[i].label)+'-'
        strr += str(L[-1].label)
        sgn = random.choice(['+','-'])
        signlist.append(sgn)
        print(strr+' ',sgn)
    return signlist

def print_cycle_pattern(cycle_list, sign_list, numbers = None):
    #print a list of cycles
    #Inputs:
    #cycle_list: list of cycles, each cycle given as list of vertices
    #sign_list: list of signs, the i-th sign belonging to the i-th cycle.
    #numbers: optional list for if one wants to show how many times we have a cycle

    if len(cycle_list)!=len(sign_list):
        raise ValueError('Length of lists is not the same')
    signlist = []
    for ii in range(len(cycle_list)):
        L = cycle_list[ii]
        strr=''
        for i in range(len(L)-1):
            strr+=str(L[i].label)+'-'
        strr += str(L[-1].label)
        sgn = sign_list[ii]
        if numbers is not None:
            print(numbers[ii],'x: '+strr+' ',sgn)
        else:
            print(strr+' ',sgn)
    return

def cycle_text(L):
    #turn a cycle into a string
    strr = ''
    for i in range(len(L) - 1):
        strr += str(L[i].label) + '-'
    strr += str(L[-1].label)
    return strr + ' '

def find_cycles(G, L, used_nodes):
    """
    subroutine that returns cycles containing the path L that do not use nodes from used_nodes
    :param G:
    :param L:
    :param used_nodes:
    :return: cycle_list
    """
    cycle_list = []
    l=set(L).union(used_nodes)
    headd = L[-1]
    taill = L[0]
    for e in headd.incidence:
        if e.tail == headd:
            if e.head == taill:
                cycle_list.append(L)
            elif e.head not in l:
                cycle_list+=find_cycles(G,L+[e.head],used_nodes)
    return cycle_list

def check_pattern(G,cycle_list, sign_list, poset=False, make_conic_basis=False, undefined_cycles=False):
    """
    #check if cycle pattern is valid
    #G- graph
    #cycle_list - list of all cycles
    #sign_list - list of signs of cycles in same order as cycle_list. Sign can be '-', '+' or 0. If undefined_cycles is true, then 0 means we don't prescribe its sign
    #poset - construct the partial order of the cycles
    #make_conic_basis - find the signed cycle basis
    #undefined_cycles - Interpret 0 as undefined instead of a 0-weight cycle
    """
    m = len(G.edges)
    E=G.edges
    V=G.vertices
    A=[]
    b=[]
    for i in range(len(cycle_list)): #construct LP
        C=cycle_list[i]
        s=sign_list[i]
        if s == '+':
            pm = -1
        elif s==0 and undefined_cycles:
            pm = 0
        else:
            pm = 1
        a = [0 for l in range(m)]
        for j in range(len(C)-1):
            for k in range(m):
                if E[k].tail == C[j] and E[k].head == C[j+1]: #if edge is in cycle
                    a[k] = pm
        for k in range(m):
            if E[k].tail == C[-1] and E[k].head == C[0]:  # if edge is last edge of cycle
                a[k] = pm
        A.append(a)
        if undefined_cycles and s==0:
            b.append(0)
        else:
            b.append(-1)
    AAB = [a for a in A]
    c = [-1*x for x in a] # so LP is not unbounded
    res = linprog(c, A, b, bounds=(None, None))
    A=AAB
    if res.status == 0:
        print('A solution was found')
        vec = res.x
        for ii in range(m):
            print(E[ii],':',vec[ii])
            E[ii].phiweight = round(vec[ii])
        if poset:
            H = poset_graph(G, A, b, c, cycle_list, sign_list)
            H.label = "Poset graph: Arrow means greater than, circle is + and square is -"
            with open('test2.gv', 'w') as filez:
                write_dot(H, filez, True)
            src3 = Source.from_file('test2.gv')
            h = Digraph()
            source_lines = str(src3).splitlines()
            source_lines.pop(0)
            source_lines.pop(-1)
            h.body += source_lines
            h.graph_attr['splines'] = 'true'
            h.graph_attr['sep'] = '1'
            h.graph_attr['scale'] = '1'
            h.render('test2.gv', format='png', cleanup=True, quiet=True, engine="neato", view=True).replace('\\', '/')
            ##Count length of conic cycle basis
            pluscycles = 0
            minuscycles = 0
            for v in H.vertices:
                if v.player == 0 and len([e for e in v.incidence if e.tail==v])==0: #+cycle without outgoing edges
                    pluscycles += 1
                elif v.player == 1 and len([e for e in v.incidence if e.head==v])==0: #-cycle without incoming edges
                    minuscycles += 1
        if make_conic_basis:
            pluscyclelist, minuscyclelist = conic_basis(G, A, b, c, cycle_list, sign_list)
            pluscycles = len(pluscyclelist)
            minuscycles = len(minuscyclelist)
            print('Basic cycles:')
            for igr in pluscyclelist + minuscyclelist:
                print(igr)
            print('G has', len(G.vertices), 'nodes and', len(G.edges), 'edges, so it has', len(G.edges)-len(G.vertices)+1, 'cycles\n in the basis if it is strongly connected')
            print('Actual number of cycles in basis:',np.linalg.matrix_rank(np.matrix(A)))
            print('Total number of cycles in graph:', len(cycle_list))
            print('Positive cycles in the conic basis:', pluscycles)
            print('Negative cycles in the conic basis:', minuscycles)
            Aneg = [A[aaaa] for aaaa in range(len(cycle_list)) if sign_list[aaaa]=='-']
            Apos = [A[aaaa] for aaaa in range(len(cycle_list)) if sign_list[aaaa] == '+']
            print('Dimension of positive cycle space:', np.linalg.matrix_rank(np.matrix(Apos)))
            print('Dimension of negative cycle space:', np.linalg.matrix_rank(np.matrix(Aneg)))
        if poset:
            return -len(H.edges)
        else:
            return -math.inf
    elif res.status == 2:
        print('cycle pattern not valid')
        # make LP smaller
        AA = [aa for aa in A] # new A
        bb = [bbb for bbb in b]
        indices = [iii for iii in range(len(b))]
        count = 0
        while count < len(indices):
            AAA = [AA[iii] for iii in range(len(indices)) if iii != count]
            bbb = [bb[iii] for iii in range(len(indices)) if iii != count]
            res2 = linprog(c, AAA, bbb, bounds=(None, None))
            if res2.status == 2:
                AA=AAA
                bb=bbb
                indices = [indices[iii] for iii in range(len(indices)) if iii != count]
            else:
                count +=1
        #print smallest set of invalid cycles
        print('Minimal set of cycles that is invalid:')
        for ind in indices:
            L=cycle_list[ind]
            strr=''
            for i in range(len(L)-1):
                strr+=str(L[i].label)+'-'
            strr += str(L[-1].label)
            sgn = sign_list[ind]
            sign_list.append(sgn)
            print(strr+' ',sgn)
        return len(indices)
    else:
        print('???????')
        return -1

def construct_big_witness_graph(N):
    """
    Construct a graph in a family that, for specific cycle patterns,
    only has exponential size witnesses
    :param N: number of segments of the right part
    :return: graph G_N, xev, xod
    xev is the vertex with the largest even label, and xod the vertex with the largest odd label
    """
    global G
    scalef = 12
    G = MyGraph(True, 0, False)
    #starting part
    x = Vertex(G, 1, False, 0, pos='0,0')
    G.add_vertex(x)
    y = Vertex(G, 2, False, 0, pos='1,1')
    G.add_vertex(y)
    x2 = Vertex(G, 3, False, 0, pos='2,0')
    G.add_vertex(x2)
    G.add_edge(Edge(x,y,-2**(N+4)-2))
    G.add_edge(Edge(y,x2,0))
    G.add_edge(Edge(x,x2,2**(N+4)-2))
    G.add_edge(Edge(x2,x,0))
    y2 = Vertex(G, 4, False, 0, pos='3,1')
    G.add_vertex(y2)
    x3 = Vertex(G, 5, False, 0, pos='4,0')
    G.add_vertex(x3)
    G.add_edge(Edge(x2, y2, -2**(N+4)))
    G.add_edge(Edge(y2, x3, 0))
    G.add_edge(Edge(x2, x3, 2**(N+4)))
    G.add_edge(Edge(x3, x, 0))
    x4 = Vertex(G, 6, False, 0, pos='6,2')
    G.add_vertex(x4)
    x5 = Vertex(G, 7, False, 0, pos='6,-2')
    G.add_vertex(x5)
    G.add_edge(Edge(x3,x4,0))
    G.add_edge(Edge(x3,x5,0))
    xprevev=x4
    xprevod=x5
    for n in range(N):
        if n==N-1 and False:
            xdoubleev = Vertex(G, 4*n+6, False, 0, pos=str(2*n+6+0.2)+',2')
            G.add_vertex(xdoubleev)
            xdoubleod = Vertex(G, 4 * n + 7, False, 0, pos=str(2*n + 6 + 0.2) + ',-2')
            G.add_vertex(xdoubleod)
            G.add_edge(Edge(xprevev,xdoubleev,0))
            xprevev = xdoubleev
            G.add_edge(Edge(xprevod,xdoubleod,0))
            xprevod=xdoubleod
        xnewev = Vertex(G, 4*n+10, False, 0, pos=str(2*n+8)+',2')
        G.add_vertex(xnewev)
        xnewod = Vertex(G, 4*n+11, False, 0, pos=str(2*n+8)+',-2')
        G.add_vertex(xnewod)
        yev = Vertex(G, 4*n+8, False, 0, pos=str(2*n+7)+',3')
        G.add_vertex(yev)
        yod = Vertex(G, 4*n+9, False, 0, pos=str(2*n+7)+',-1')
        G.add_vertex(yod)
        G.add_edge(Edge(xprevev,yev,-2**(N+5+n)+2**(n+2)))
        G.add_edge(Edge(yev, xnewev, 0))
        G.add_edge(Edge(xprevev, xnewev, 2**(N+5+n)-2**(n+2)))
        G.add_edge(Edge(xnewev,x,0))
        G.add_edge(Edge(yod, xnewod, 0))
        if n<N-1:
            G.add_edge(Edge(xprevod, xnewod, 2**(N+5+n)+2**(n+2)))
            G.add_edge(Edge(xprevod, yod, -2 ** (N + 5 + n) - 2 ** (n + 2)))
        else:
            G.add_edge(Edge(xprevod, xnewod, 2**(N+5+n)+2**(n+2)-1)) #make odd to detect it
            G.add_edge(Edge(xprevod, yod, -2 ** (N + 5 + n) - 2 ** (n + 2)-1))
            G.add_edge(Edge(xnewod, xprevev, 0))
            G.add_edge(Edge(xnewev, xprevod, 0))
        G.add_edge(Edge(xnewod, x, 0))
        xprevev = xnewev
        xprevod = xnewod
    return G, xnewev, xnewod

def find_smallest_witness(N):
    """
    For the specific graph G_N constructed by construct_big_witness_graph(N),
    find the smallest possible witness
    :param N:
    :return: G, cycle_list, sign_list, res
    res is a scipy object containing the results and info of solving the ILP
    """
    G, endvert1, endvert2 = construct_big_witness_graph(N)
    with open('test2.gv', 'w') as filez:
        write_dot(G, filez, True)
    src3 = Source.from_file('test2.gv')
    h = Digraph()
    source_lines = str(src3).splitlines()
    source_lines.pop(0)
    source_lines.pop(-1)
    h.body += source_lines
    h.graph_attr['splines'] = 'true'
    h.graph_attr['sep'] = '1'
    h.graph_attr['scale'] = '1'
    h.render('test2.gv', format='png', cleanup=True, quiet=True, engine="neato", view=True).replace('\\', '/')

    cycle_list = make_cycle_list(G)
    sign_list = []
    for cyc in cycle_list:
        cyc_weight = 0
        for i in range(len(cyc)):
            edgeset = G.find_edge(cyc[i-1], cyc[i])
            for e in edgeset:
                if e.head == cyc[i]:
                    cyc_weight += e.weight
        if cyc_weight%2==1 and abs(cyc_weight)<2**(N+3) and endvert1 in cyc and endvert2 in cyc:
            if len(cyc) <= 6:
                if cyc_weight == 2**(N+2)-1:
                    sign_list.append('-')
                elif cyc_weight == -2**(N+2)-1:
                    sign_list.append('+')
                else:
                    raise ValueError('Wrong weight for cycle', cyc)
            else: #cycle weight depends on whether it has the 2**(N+2) or the -2**(N+2)
                if cyc_weight > 0:
                    fake_cyc_weight = cyc_weight-(2**(N+2)-1)
                else:
                    fake_cyc_weight = cyc_weight+(2**(N+2)+1)
                if fake_cyc_weight > 0:
                    sign_list.append('+')
                else:
                    sign_list.append('-')
        elif cyc_weight > 0:
            sign_list.append('+') #the one cycle gets the wrong sign
        elif cyc_weight == 0:
            sign_list.append('0') #should not really make a difference - but makes it a lot easier to deal with the ILP
            print('Warning: 0 weight cycle. ILP will break.')
        elif cyc_weight <0:
            sign_list.append('-')
    ######## write ILP:
    m = len(G.edges)
    E = G.edges
    V = G.vertices
    A = []
    b = []
    for i in range(len(cycle_list)):  # construct LP
        C = cycle_list[i]
        s = sign_list[i]
        if s == '+':
            pm = -1
        elif s == 0:
            pm = 0
        else:
            pm = 1
        a = [0 for l in range(m)]
        for j in range(len(C) - 1):
            for k in range(m):
                if E[k].tail == C[j] and E[k].head == C[j + 1]:  # if edge is in cycle
                    a[k] = pm
        for k in range(m):
            if E[k].tail == C[-1] and E[k].head == C[0]:  # if edge is last edge of cycle
                a[k] = pm
        A.append(a)
        if s == 0:
            b.append(0)
        else:
            b.append(-1)
    A=np.array(A)
    A=np.transpose(A)
    b=np.array(b)
    print('Generated graph with',len(G.vertices),'nodes and',len(G.edges),'edges')
    print('Number of cycles:', len(cycle_list))
    print('Size of ILP matrix:',np.shape(A))

    #A is multiplied by -1 here, but it does not matter
    c = np.array([1 for i in range(len(cycle_list))])
    Ab = np.array([0 for i in range(m)])
    bb = np.array([-1])
    const1 = LinearConstraint(A, Ab, Ab)
    const2 = LinearConstraint(b, ub=bb)

    res = milp(c, integrality=c, constraints=[const1,const2])
    print('Scipy ILP optimization status:',res.status)
    print('Optimal value:', res.fun)
    print('Optimum:', res.x)
    optimum = list(res.x)
    critical_cycles = []
    critical_signs = []
    numbers = []
    for i in range(len(optimum)):
        if optimum[i]>0.5:
            critical_cycles.append(cycle_list[i])
            critical_signs.append(sign_list[i])
            numbers.append(round(optimum[i]))
    print('Minimal witness:')
    print_cycle_pattern(critical_cycles,critical_signs, numbers)

    return G, cycle_list, sign_list, res

def conic_basis(G, A, b, c, cycle_list, sign_list):
    """
    Subroutine for finding the smallest +-cycles and the smallest --cycles
    :param G:
    :param A:
    :param b:
    :param c:
    :param cycle_list:
    :param sign_list:
    :return: pluscycs, mincycs : The smallest +-cycles resp. biggest --cycles
    """
    N = len(cycle_list)
    H = MyGraph(True, N, True)
    for ii in range(N):
        L = cycle_list[ii]
        strr = ''
        for i in range(len(L) - 1):
            strr += str(L[i].label) + '-'
        strr += str(L[-1].label)
        strr += ' '+sign_list[ii]
        H.vertices[ii].label = strr
        H.vertices[ii].player = (sign_list[ii] != '+')
    AA = [aa for aa in A]
    bb = [bbb for bbb in b]
    pluscycs = [] #positive basic cycles
    mincycs = [] #negative basic cycles
    for cyc in range(N):
        AAA = [i for i in AA]
        AAA[cyc]=[-AAA[cyc][jk] for jk in range(len(AAA[cyc]))]
        res2 = linprog(c, AAA, bb, bounds=(None, None))
        if res2.status == 0:  # when startt > endd
            if sign_list[cyc]=='+':
                pluscycs.append(H.vertices[cyc].label)
            else:
                mincycs.append(H.vertices[cyc].label)
    return pluscycs, mincycs

def poset_graph(G, A, b, c, cycle_list, sign_list):
    """
    Subroutine to construct the partially ordered set of the cycles that follows from the cycle pattern
    :param G: input graph
    :param A: matrix related to realizability cone
    :param b: right hand side related to realizability cone
    :param c: objective function
    :param cycle_list: list of cycles
    :param sign_list: list of signs
    :return: Partially ordered set of cycles given as a graph: edge (C_1,C_2) means that C_1 is bigger than C_2
    """
    N = len(cycle_list)
    H = MyGraph(True, N, True)
    for ii in range(N):
        L = cycle_list[ii]
        strr = ''
        for i in range(len(L) - 1):
            strr += str(L[i].label) + '-'
        strr += str(L[-1].label)
        H.vertices[ii].label = strr
        H.vertices[ii].player = (sign_list[ii] != '+')
    AA = [aa for aa in A]+[[0 for i in range(len(A[0]))]]
    bb = [bbb for bbb in b]+[-1]
    for startt in range(N):
        for endd in range(N):
            if startt != endd: #test if cycle start is larger than cycle end
                if not (sign_list[startt] == '+' and sign_list[endd] == '-'):
                    extraeq = [abs(AA[startt][i])-abs(AA[endd][i]) for i in range(len(AA[startt]))]
                    AA[-1] = extraeq
                    res2 = linprog(c, AA, bb, bounds=(None, None))
                    if res2.status == 2: #when startt > endd
                        H.add_edge(Edge(H.vertices[startt],H.vertices[endd]))
    for v in H.vertices:
        reached = set() #list of vertices reachable via a k-long path with k>=2
        queue = [] #queue of to-check vertices
        for e in v.incidence:
            if e.tail == v:
                queue.append(e.head)
        while queue != []:
            q = queue[0]
            for e in q.incidence:
                if e.tail == q and e.head not in queue and e.head not in reached:
                    queue.append(e.head)
                    reached.add(e.head)
            queue.pop(0)
        for e in v.incidence:
            if e.tail == v and e.head in reached:
                H.remove_edge(e)
    return H

def check_parity_valid(G,L, SL):
    """
    check if there is a parity game that has cycles L and sign pattern SL
    """
    E = G.edges
    m=len(E)
    prios = [0 for i in range(len(E))]
    EL = [[] for j in range(len(E))] #list of cycle numbers per edge
    for i in range(len(L)):  # construct LP
        C = L[i]
        for j in range(len(C) - 1):
            for k in range(m):
                if E[k].tail == C[j] and E[k].head == C[j + 1]:  # if edge is in cycle
                    EL[k].append(i)
        for k in range(m):
            if E[k].tail == C[-1] and E[k].head == C[0]:  # if edge is last edge of cycle
                EL[k].append(i)
    remaining_cycles = set(range(len(L))) #cycles that still have no assigned edges
    current_prio = 2*len(E)
    remaining_edges = set(range(len(E)))
    while len(remaining_cycles)>0:
        changed = False
        for e in remaining_edges:
            pluspos = True #if edge can be large even
            minpos = True #if edge can be large odd
            for f in EL[e]:
                if f in remaining_cycles and SL[f]=='+':
                    minpos = False
                elif f in remaining_cycles and SL[f]=='-':
                    pluspos = False
            if pluspos:
                prios[e]=current_prio
                current_prio -= 2
                changed = True
                remaining_edges.remove(e)
                remaining_cycles = remaining_cycles.difference(EL[e])
                break
            elif minpos:
                prios[e]=current_prio-1
                current_prio -= 2
                changed = True
                remaining_edges.remove(e)
                remaining_cycles = remaining_cycles.difference(EL[e])
                break
        if not changed:
            return False, 0
    return True, prios


if __name__=='__main__':
    ###Uncomment to generate random graph and check if its cycle pattern is realizable
    #random.seed(1)
    #np.random.seed(1)
    #G = make_random_graph(8,18)
    #L=make_cycle_list(G)
    #SL = generate_and_print_random_sign_list(L) ##SL=generate_and_print_random_valid_sign_list(L) to guarantee realizability
    #sizez = check_pattern(G, L, SL)
    ###Set to true to save image of graph
    if False:
        G.label = "Graph"
        with open('test.gv', 'w') as filez:
            write_dot(G, filez, True)
        src2 = Source.from_file('test.gv')
        g = Digraph()
        source_lines = str(src2).splitlines()
        source_lines.pop(0)
        source_lines.pop(-1)
        g.body += source_lines
        g.graph_attr['splines'] = 'true'
        g.graph_attr['sep'] = '1'
        g.graph_attr['scale'] = '1'
        g.render('test.gv', format='png', cleanup=True, quiet=True, engine="dot", view=True).replace('\\', '/')

    G, cycle_list, sign_list, res = find_smallest_witness(10)
    #print('Cycle pattern:')
    #print_cycle_pattern(cycle_list, sign_list)
    #check_pattern(G, cycle_list, sign_list, 0, poset=False)

    #help(check_parity_valid)