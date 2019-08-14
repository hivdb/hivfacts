#! /usr/bin/env python
from hivfacts import HIVAPOBEC


def main():
    apobec = HIVAPOBEC()
    sql = [
        'CREATE TABLE IF NOT EXISTS `tblLUAPOBEC2` (',
        "  `Gene` enum('PR','RT','IN') NOT NULL,",
        '  `Position` smallint(3) unsigned NOT NULL,',
        '  `AA` char(1) NOT NULL,',
        "  `DRM` enum('Yes','No') NOT NULL,",
        '  PRIMARY KEY (`Gene`,`Position`,`AA`)',
        ') ENGINE=MyISAM DEFAULT CHARSET=utf8;',
        '',
        'TRUNCATE TABLE `tblLUAPOBEC2`;'
        '',
        '',
    ]

    for mut in apobec.get_apobec_list():
        sql.append(
            "INSERT INTO `tblLUAPOBEC2` (Gene, Position, AA, DRM) VALUES "
            "('{gene}', {position}, '{aa}', 'No');"
            .format(**mut))

    for mut in apobec.get_apobec_drm_list():
        sql.append(
            "INSERT INTO `tblLUAPOBEC2` (Gene, Position, AA, DRM) VALUES "
            "('{gene}', {position}, '{aa}', 'Yes');"
            .format(**mut))

    print('\n'.join(sql))


if __name__ == '__main__':
    main()
