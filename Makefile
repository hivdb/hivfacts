refresh-yamls:
	@pipenv run python scripts/yaml2json.py

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
