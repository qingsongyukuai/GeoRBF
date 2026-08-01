FROM python:3.12.3-slim-bookworm@sha256:afc139a0a640942491ec481ad8dda10f2c5b753f5c969393b12480155fe15a63

WORKDIR /oracle

COPY requirements.lock ./
RUN python -m pip install --require-hashes -r requirements.lock

COPY . ./

CMD ["sh", "-ec", "python -m unittest discover -s tests -v && python generate.py --check && python generate.py --check && python verify.py"]
