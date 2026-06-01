import http.server, ssl, socket

port = 8443
handler = http.server.SimpleHTTPRequestHandler
httpd = http.server.HTTPServer(('0.0.0.0', port), handler)
ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
ctx.load_cert_chain('cert.pem', 'key.pem')
httpd.socket = ctx.wrap_socket(httpd.socket, server_side=True)

ip = socket.gethostbyname(socket.gethostname())
print(f"Serving on https://{ip}:{port}")
print("Open that URL on your phone — accept the certificate warning")
httpd.serve_forever()
