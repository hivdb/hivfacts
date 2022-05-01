data/drms_%.json: data/drms_%.csv
	@pipenv run python scripts/drms_csv2json.py $<

data/mutation-type-pairs_%.json: data/mutation-type-pairs_%.csv
	@pipenv run python scripts/csv2json.py $<

data/apobec%/apobecs.json: data/apobec%/apobecs.csv
	@pipenv run python scripts/csv2json.py $<

data/apobec%/apobec_drms.json: data/apobec%/apobec_drms.csv
	@pipenv run python scripts/csv2json.py $<

data/aapcnt/rx-all_subtype-%.json: data/aapcnt/rx-all_subtype-%.csv
	@pipenv run python scripts/csv2json.py $<

data/%.json: data/%.yml
	@pipenv run python scripts/yaml2json.py $<

data: data/*.json data/apobecs/*.json data/apobecs-hiv2/*.json

refresh-yamls:
	@pipenv run python scripts/yaml2json.py

refresh-csvs:
	@scripts/update-csvs.sh

refresh-drmlist:
	@scripts/update-drmlist.sh

refresh-sdrmlist:
	@scripts/update-sdrmlist.sh

refresh-tsmlist:
	@scripts/update-tsmlist.sh

refresh-muttypepairs:
	@scripts/update-muttypepairs.sh

refresh-aapcnt:
	@scripts/update-aapcnt.sh

refresh-codonpcnt:
	@scripts/update-codonpcnt.sh

refresh-condcomments:
	@scripts/update-condcomments.sh

refresh-apobec:
	@scripts/prepare-apobec-data.sh
	@scripts/update-apobec.sh

refresh-data: refresh-drmlist refresh-sdrmlist refresh-tsmlist refresh-muttypepairs

build-java:
	@cd hivfacts-java; ./gradlew build

release-java:
	@cd hivfacts-java; ./gradlew build bintrayUpload

build-python:
	@cd hivfacts-python; ./setup.py copy_data sdist bdist_wheel

release-python:
	@cd hivfacts-python; ./setup.py copy_data sdist bdist_wheel upload

.PHONY: release-java
