# BiochemistBuddy Primer3 Engine

An isolated GPLv2 service that exposes Primer3 and `primer3-py` to separately
licensed applications over an authenticated HTTP boundary.

The service is intended to run inside the same customer-controlled VPC as its
caller. It does not persist inputs or results, call external services, or include
BiochemistBuddy product code. Customer sequences remain inside the customer silo.

## API

- `GET /healthz`
- `POST /v1/design/pcr`

The design endpoint requires a short-lived BiochemistBuddy service JWT with:

- `customer_id` matching `DEPLOYMENT_CUSTOMER_ID`
- scope `primer-engine:execute`
- issuer `biochemistbuddy`
- audience `biochemistbuddy-platform`

`SERVICE_JWT_SECRET` must be shared securely with the calling service inside the
customer deployment.

## Development

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
pytest
ruff check .
docker build -t biochemistbuddy-primer3-engine .
```

All automated tests use synthetic DNA fixtures.

## License

This repository is licensed under GPLv2 only. It includes and links to
GPLv2-licensed Primer3 through `primer3-py`. See `LICENSE` and retain the source,
copyright notices, build files, and license when distributing an image.
