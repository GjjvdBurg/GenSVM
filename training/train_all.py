#!/usr/bin/env python

import os

def main():
    for f in os.listdir('.'):
        if f.endswith('.training'):
            os.system('../trainMSVMMajdataset %s > %s.output' % (f, f))

if __name__ == '__main__':
    main()
