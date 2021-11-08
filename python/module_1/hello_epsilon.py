from numpy import finfo

if __name__ == "__main__":
    # According to https://scipy-lectures.org/_downloads/ScipyLectures-simple.pdf
    # python has int, float, complex and booleans
    print("Machine epsilon: %s" % finfo(float).eps)
