import agrparse
from math import sqrt, asin, sin, pi

def lamberts_equation(a, s, c):

    alph = 2*pi - 2*asin(sqrt(s/(2*a)))
    beta = 2*asin(sqrt((s-c)/(2*a)))
    tof = a**(3/2) * (alph - beta - (sin(alph)-sin(beta)))

    return tof

def main():

    parser = argparse.ArgumentParser(description="A simple Python script that shows lamberts problem trajectories")
    
    # Add command line arguments
    parser.add_argument('--example', 
                        type = str, 
                        choices = ['v', 'm', 'j'],
                        help = "Which Lambert's example you want to use, out of Mars (m), Jupiter (j) or Venus (v) transfer")
    
    # parser.add_argument('--output', type=str, help="Output file")
    # parser.add_argument('--verbose', action='store_true', help="Enable verbose mode")
    
    args = parser.parse_args()
    
    # Access the argument values
    example = args.example