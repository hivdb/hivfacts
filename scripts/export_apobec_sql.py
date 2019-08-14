#! /usr/bin/env python
from hivfacts import HIVAPOBEC


def main():
    apobec = HIVAPOBEC()
    sql = [
        'CREATE TABLE IF NOT EXISTS `tblLUAPOBEC` (',
        "  `Gene` enum('PR','RT','IN') NOT NULL,",
        '  `Position` smallint(3) unsigned NOT NULL,',
        '  `AA` char(1) NOT NULL,',
        '  PRIMARY KEY (`Gene`,`Position`,`AA`)',
        ') ENGINE=MyISAM DEFAULT CHARSET=utf8;',
        '',
        'TRUNCATE TABLE `tblLUAPOBEC`;'
        '',
        '',
    ]

    for mut in apobec.get_apobec_list():
        sql.append(
            "INSERT INTO `tblLUAPOBEC` (Gene, Position, AA) VALUES "
            "('{gene}', {position}, '{aa}');"
            .format(**mut))

    print('\n'.join(sql))


if __name__ == '__main__':
    main()
