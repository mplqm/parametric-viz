import http.server, ssl, os

PORT     = 8444
USDZ_OUT = 'protein.usdz'

class Handler(http.server.BaseHTTPRequestHandler):

    def do_OPTIONS(self):
        self.send_response(200)
        self._cors()
        self.end_headers()

    def do_POST(self):
        if self.path != '/save-usdz':
            self.send_response(404); self.end_headers(); return
        length = int(self.headers.get('Content-Length', 0))
        data   = self.rfile.read(length)
        with open(USDZ_OUT, 'wb') as f:
            f.write(data)
        self.send_response(200)
        self._cors()
        self.send_header('Content-Type', 'application/json')
        self.end_headers()
        self.wfile.write(b'{"ok":true}')
        print(f'  saved {len(data)} bytes → {USDZ_OUT}')

    def do_GET(self):
        if self.path != '/protein.usdz':
            self.send_response(404); self.end_headers(); return
        if not os.path.exists(USDZ_OUT):
            self.send_response(404); self.end_headers(); return
        data = open(USDZ_OUT, 'rb').read()
        self.send_response(200)
        self._cors()
        self.send_header('Content-Type',   'model/vnd.usdz+zip')
        self.send_header('Content-Length', str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _cors(self):
        self.send_header('Access-Control-Allow-Origin',  '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')

    def log_message(self, fmt, *args):
        print(f'  {self.address_string()}  {fmt % args}')

httpd = http.server.HTTPServer(('0.0.0.0', PORT), Handler)
ctx   = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
ctx.load_cert_chain('cert.pem', 'key.pem')
httpd.socket = ctx.wrap_socket(httpd.socket, server_side=True)
print(f'USDZ server → https://0.0.0.0:{PORT}')
print(f'  POST /save-usdz     save body as {USDZ_OUT}')
print(f'  GET  /protein.usdz  serve with model/vnd.usdz+zip')
httpd.serve_forever()
