FROM public.ecr.aws/lambda/python:3.12

LABEL org.opencontainers.image.source="https://github.com/LibreTexts/molview-lambda-api"

COPY requirements.txt ${LAMBDA_TASK_ROOT}

RUN pip install -r requirements.txt

COPY api.py ${LAMBDA_TASK_ROOT}

CMD [ "api.handler" ]