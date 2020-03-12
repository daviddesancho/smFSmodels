import numpy as np
def calc_life(trajs, ub=5, lb=-5):
    """
    Identifies transition paths and returns lifetimes of states.
    
    Parameters
    ----------
    trajs : list of lists
        Set of trajectories.
        
    ub, lb : float
        Cutoff value for upper and lower states.
    """
    try:
        assert ub > lb
    except AssertionError:
        print (" Upper bound is lower than lower bound")
        return
    
    lifeA = []
    lifeB = []
    time = 0
    for tr in trajs:
        state = None
        ntp = 0
        time_prev = 0
        for t,q in enumerate(tr):
            # assign state when beyond boundaries
            if q > ub: # state "B"
                if state == 'A':
                    ntp +=1
                    lifeA.append(time - time_prev)
                    time_prev = time
                state = 'B'
            elif q < lb: # state "A"
                if state == 'B':
                    ntp +=1
                    lifeB.append(time - time_prev)
                    time_prev = time
                state = 'A'
            else:
                if state == 'A' and q < ub:
                    time = t
                elif state == 'B' and q > lb:
                    time = t
    return lifeA, lifeB 


def calc_life_multi(trajs, bounds=[[-3,-1], [1,3], [6,8]]):
    """
    Identifies transition paths and returns lifetimes of states.
    
    Parameters
    ----------
    trajs : list of lists
        Set of trajectories.
        
    bounds : list 
        Limits for states.
    """
    
    life = [[],[],[]]
    tau = {} 
    time = 0
    for tr in trajs:
        state = None
        ntp = 0
        for t,q in enumerate(tr):
            # assign state when beyond boundaries
            for i,b in enumerate(bounds):
                if b[0] < q < b[1]:
                    if state != i:
                        try:
                            life[state].append(time - time_prev)
                            tau[i, state].append(time - time_prev)
                        except TypeError:
                            pass
                        except KeyError:
                            tau[(i, state)] = [time - time_prev]
                        state = i
                        time_prev = t
                        break
                    state = i
                time = t

    return life, tau

def lifetimes(data, f_bound=0, u_bound=10):
    """
    Estimates lifetimes using a transition path analysis. Transitions are
    only assigned from one state to the other when the core of the other
    state is reached.
    
    Parameters
    ----------
    data : np.array
        Time series data containing times and corrected extensions.
        
    Returns
    -------
    tau_f : list
        Waiting times in the unfolded state.
        
    tau_u : list
        Waiting times in the folded state.        
        
    data_f : list
        Stretches of data corresponding to the folded segments.
        
    data_u : list
        Stretches of data corresponding to the unfolded segments.
    
    tp_f : list
        Transition paths for folding.
        
    tp_u : list
        Transition paths for unfolding.  
    """
    folded = False
    unfolded = False
    maybetp = []
    recrossings = []
    t = data[:,0]
    dist = data[:,1]
    data_f = []
    data_u = []
    tau_u = []
    tau_f = []
    tp_f = []
    tp_u = []
    time = 0
    data_cum = []
    for t, d in zip(t,dist):
        if d <= f_bound:
            folded = True
            if unfolded:
                #print ' Refolding event: %g'%t,
                tp_u.append(np.array(maybetp))
                tau_f.append(t-time)
                data_u.append(np.array(data_cum))
                unfolded = False
                time = t
                data_cum = []
            else:
                recrossings.append(np.array(maybetp))
                for tt,dd in maybetp:
                    data_cum.append([tt,dd])
            data_cum.append([t,d])

            maybetp = []

        elif u_bound <= d :
            unfolded = True
            if folded:
                tp_f.append(np.array(maybetp))
                #print ' Unfolding event: %g'%t,
                tau_u.append(t-time)
                data_f.append(np.array(data_cum))
                folded = False
                time = t
                data_cum = []
            else:
                recrossings.append(np.array(maybetp))
                for tt,dd in maybetp:
                    data_cum.append([tt,dd])
            data_cum.append([t,d])
            maybetp = []
        else:
            maybetp.append([t,d])
        
    if unfolded:
        #print ' Refolding event: %g'%t,
        tau_f.append(t-time)
        data_u.append(np.array(data_cum))

    if folded:
#        tp_f.append(np.array(maybetp))
        #print ' Unfolding event: %g'%t,
        tau_u.append(t-time)
        data_f.append(np.array(data_cum))
    return tau_f, tau_u, data_f, data_u, tp_f, tp_u, recrossings
