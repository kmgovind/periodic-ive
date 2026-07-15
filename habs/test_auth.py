import earthaccess

# This will prompt you to log in once.
# It saves your credentials to a .netrc file so you don't have to do this again.
auth = earthaccess.login(strategy="interactive")

if auth.authenticated:
    print("Success! You're ready to haunt NASA's servers.")
else:
    print("Authentication failed. Double-check your Earthdata username/password.")