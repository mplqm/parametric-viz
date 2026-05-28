"""
HTTPS server for local WebXR AR testing.
WebXR requires a secure context (HTTPS or localhost).

Usage:
    python3 serve_https.py          # generates cert on first run
    python3 serve_https.py 8443     # optional port argument

Testing on iPhone (same local network):
    1. Run this script — note the LAN IP printed
    2. Open https://<LAN_IP>:8443 in Safari on iPhone
    3. Accept the security warning (tap Advanced → Proceed)
    4. Settings → Safari → Advanced → Feature Flags → WebXR (enable if needed)
    5. Tap the AR button — point at a table, tap to place

Desktop fallback:
    http://localhost:8080  (standard HTTP, no AR button)
    https://localhost:8443 (HTTPS, AR button shows if device supports it)
"""

import http.server
import ssl
import os
import socket
import subprocess
import sys

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8443
CERT = "cert.pem"
KEY  = "key.pem"

# Generate self-signed cert if not already present
if not (os.path.exists(CERT) and os.path.exists(KEY)):
    print("Generating self-signed TLS certificate…")
    subprocess.run([
        "openssl", "req", "-x509", "-newkey", "rsa:2048",
        "-keyout", KEY, "-out", CERT,
        "-days", "365", "-nodes",
        "-subj", "/CN=localhost",
    ], check=True)
    print(f"  {CERT} + {KEY} written\n")

# Resolve LAN IP for mobile access instructions
try:
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    lan_ip = s.getsockname()[0]
    s.close()
except Exception:
    lan_ip = "127.0.0.1"

print(f"Serving HTTPS on port {PORT}")
print(f"  Desktop : https://localhost:{PORT}")
print(f"  Mobile  : https://{lan_ip}:{PORT}")
print(f"\nNote: your browser will show a security warning for the self-signed")
print(f"certificate — click 'Advanced' → 'Proceed to site' to continue.\n")

httpd   = http.server.HTTPServer(("", PORT), http.server.SimpleHTTPRequestHandler)
context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
context.load_cert_chain(CERT, KEY)
httpd.socket = context.wrap_socket(httpd.socket, server_side=True)

try:
    httpd.serve_forever()
except KeyboardInterrupt:
    print("\nServer stopped.")
