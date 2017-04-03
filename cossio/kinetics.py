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
        print " Upper bound is lower than lower bound"
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
