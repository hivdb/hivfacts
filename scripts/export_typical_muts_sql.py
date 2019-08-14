#! /usr/bin/env python
from hivfacts import HIVAAPcnt
from hivdbql.hiv_references import HIV_REFERENCES


def main():
    aapcnt = HIVAAPcnt('all', 'all')
    sql = [
        'CREATE TABLE IF NOT EXISTS `tblLUTypicalMutations` (',
        "  `Gene` enum('PR','RT','IN') NOT NULL,",
        '  `Position` smallint(5) unsigned NOT NULL,',
        '  `Cons` char(1) NOT NULL,',
        '  `AA` char(1) NOT NULL,',
        '  `Perc` float NOT NULL,',
        '  `Total` mediumint(9) NOT NULL,',
        '  PRIMARY KEY (`Gene`,`Position`,`AA`),',
        '  KEY `Gene` (`Gene`,`Position`)',
        ') ENGINE=MyISAM DEFAULT CHARSET=utf8;',
        '',
        'TRUNCATE TABLE `tblLUTypicalMutations`;',
        '',
        '',
    ]

    for mut in aapcnt.get():
        if mut['isUnusual']:
            continue
        cons = (
            HIV_REFERENCES['CON_{}'.format(mut['gene'])]
            [0][mut['position'] - 1]
        )
        if mut['aa'] == cons:
            continue
        sql.append(
            'INSERT INTO `tblLUTypicalMutations` '
            "(Gene, Position, Cons, AA, Perc, Total) VALUES "
            "('{gene}', {position}, '{cons}', '{aa}', {pcnt}, {total});"
            .format(**mut, cons=cons, pcnt=mut['percent'] * 100))

    print('\n'.join(sql))


if __name__ == '__main__':
    main()
