FROM python:3.12-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

RUN mkdir -p /app/tmp

EXPOSE 5590

CMD ["gunicorn", "-b", "0.0.0.0:5590", "-w", "2", "--timeout", "120", "app.main:app"]
