#! /usr/bin/env python

from hivfacts import HIVAAPcnt, HIVAPOBEC


def main():
    apobecs = HIVAPOBEC()
    aapcnts = HIVAAPcnt('all', 'all')
    matrix = [[0, 0], [0, 0]]

    for mutation in aapcnts:
        is_apobec = apobecs.is_apobec_mutation(
            mutation['gene'],
            mutation['position'],
            mutation['aa']
        )
        is_unusual = mutation['isUnusual']
        matrix[is_apobec][is_unusual] += 1

    print('               (+) APOBEC    (-) APOBEC')
    print('(+) Unusual   ',
          '{: 11}'.format(matrix[1][1]),
          '    ',
          '{: 10}'.format(matrix[0][1]),
          sep='')
    print('(-) Unusual   ',
          '{: 11}'.format(matrix[1][0]),
          '    ',
          '{: 10}'.format(matrix[0][0]),
          sep='')


if __name__ == '__main__':
    main()
