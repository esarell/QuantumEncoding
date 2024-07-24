import ProbabilityDistribution as PD
import GroverRudolph_LPW.GroverRudolph as GR

def TestSimpleDistribution():
    '''
    Testing the GR algorithm for a straight line
    Seeing if this can be encoded into the amplitude
    '''
    distribution = PD.straightLine(3)
    print(distribution)
    amplitude = GR.Grover_Rudolph_func_small(3,distribution)
    print(amplitude)
    if amplitude == distribution:
        return True
    else:
        return False


if __name__ == "__main__":
    print("Simple Test:"+str(TestSimpleDistribution()))
    