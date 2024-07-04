# MolView (Lambda) API
The code in this repo is used to build a Docker image (based on the AWS Lambda Python base image) that exports a Lambda
event handler to compare or validate diagrams made in MolView.

## Building the Image
```shell
docker build --platform linux/amd64 -t molview-lambda-api:VERSION . 
```

## Usage
The Docker image is designed to be used as a Lambda function, particularly behind API Gateway. The following paths
should be mapped to the Lambda function integration:

- `/api/v1/compare`
- `/api/v1/validate`

The API requires providing an `API_KEY` environment variable that must match the content of a Bearer token in the
`Authorization` header. The API always returns a JSON body.