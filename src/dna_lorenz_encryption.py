"""
DNA-based image encryption using MD5 and the Lorenz system (DNA_LorenzMD5).

University of Khartoum - Faculty of Mathematical Sciences and Informatics
B.Sc. graduation project, February 2023.

Usage:
    python dna_lorenz_encryption.py
    python dna_lorenz_encryption.py --image 77.png --key mysecret [--plot]

Outputs (saved in the current directory):
    encrypt_<name>.png    encrypted image
    recovered_<name>.png  decrypted image
"""

import argparse
import hashlib
import os
import textwrap
from bisect import bisect_left as bsearch

import cv2
import numpy as np
from PIL import Image
from scipy.integrate import odeint

# ---------------------------------------------------------------------------
# Global constants
# ---------------------------------------------------------------------------

# Lorenz parameters: sigma = 10, beta = 8/3, rho = 28 (chaotic regime)
LORENZ_A, LORENZ_B, LORENZ_C = 10, 8 / 3, 28

# Maximum time point of the Lorenz integration
TMAX = 100

# DNA encoding rule: A = 00, T = 01, G = 10, C = 11
dna = {}
dna["00"] = "A"
dna["01"] = "T"
dna["10"] = "G"
dna["11"] = "C"

# DNA decoding
dna["A"] = [0, 0]
dna["T"] = [0, 1]
dna["G"] = [1, 0]
dna["C"] = [1, 1]

# DNA XOR
dna["AA"] = dna["TT"] = dna["GG"] = dna["CC"] = "A"  # 00 xor 00 => 00
dna["AG"] = dna["GA"] = dna["TC"] = dna["CT"] = "G"  # 01 xor 11 => 10
dna["AC"] = dna["CA"] = dna["GT"] = dna["TG"] = "C"  # 00 xor 11 => 11
dna["AT"] = dna["TA"] = dna["CG"] = dna["GC"] = "T"  # 10 xor 11 => 01


# ---------------------------------------------------------------------------
# Lorenz system
# ---------------------------------------------------------------------------

def lorenz(X, t, a, b, c):
    """The three Lorenz equations: dx/dt, dy/dt, dz/dt."""
    x, y, z = X
    x_dot = -a * (x - y)
    y_dot = c * x - y - x * z
    z_dot = -b * z + x * y
    return x_dot, y_dot, z_dot


# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------

IMAGE_EXTENSIONS = (".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff", ".webp")


def find_image(name):
    """
    Return the path of an image. The extension is optional:
    '77' finds '77.png' (or .jpg, .jpeg, ...).
    """
    name = name.strip().strip('"').strip("'")
    if not name:
        raise ValueError("Error: no image name given.")
    if os.path.isfile(name):
        return name

    exts = IMAGE_EXTENSIONS + tuple(e.upper() for e in IMAGE_EXTENSIONS)
    matches = [name + e for e in exts if os.path.isfile(name + e)]
    if not matches:
        raise FileNotFoundError(
            "Image not found: '{0}' (current folder: {1})".format(name, os.getcwd())
        )
    if len(matches) > 1:
        print("Several matches found {0}, using {1}".format(matches, matches[0]))
    return matches[0]


def image_input_from_user():
    """Ask the user for the image name (the extension is optional)."""
    path = find_image(input("Please enter file name:  "))
    print("Image loaded: " + path)
    return path


def key_input_from_user():
    """Ask the user for the secret key."""
    return input("Enter key: ")


def split_into_rgb_channels(image):
    """Split an OpenCV (BGR) image into its R, G, B arrays."""
    red = image[:, :, 2]
    green = image[:, :, 1]
    blue = image[:, :, 0]
    return red, green, blue


# ---------------------------------------------------------------------------
# Key generation
# ---------------------------------------------------------------------------

def securekey(file_path, key):
    """
    Build the MD5 key from the user secret and the image pixels.
    Returns (hex digest, width, height).
    """
    img = Image.open(file_path).convert("RGB")
    m, n = img.size  # m = width, n = height
    print("pixels: {0}  width: {1} height: {2}".format(m * n, m, n))
    pix = img.load()

    plainimage = []
    for y in range(n):
        for x in range(m):
            for k in range(0, 3):
                plainimage.append(pix[x, y][k])

    key_in_bytes = key.encode()
    digest = hashlib.md5(key_in_bytes)       # 128-bit MD5 (32 hex characters)
    digest.update(bytearray(plainimage))     # feed the image data
    return digest.hexdigest(), m, n


def update_lorenz(key):
    """
    Derive the Lorenz initial conditions (x0, y0, z0) from the hash.
    Pure function: the same key always gives the same result.
    """
    # hex -> integer -> 256-bit binary string (zero padded on the left)
    key_bin = bin(int(key, 16))[2:].zfill(256)

    # split into 32 parts of 8 bits: k1 ... k32
    parts = textwrap.wrap(key_bin, 8)
    k = {"k{0}".format(i + 1): p for i, p in enumerate(parts)}

    t1 = t2 = t3 = 0
    for i in range(1, 12):
        t1 ^= int(k["k{0}".format(i)], 2)
    for i in range(12, 23):
        t2 ^= int(k["k{0}".format(i)], 2)
    for i in range(23, 33):
        t3 ^= int(k["k{0}".format(i)], 2)

    return t1 / 256, t2 / 256, t3 / 256


# ---------------------------------------------------------------------------
# DNA coding
# ---------------------------------------------------------------------------

def decompose_matrix(file_path):
    """Read the image and return its B, G, R channels as matrices."""
    image = cv2.imread(file_path)
    if image is None:
        raise FileNotFoundError("OpenCV could not read the image: " + file_path)
    red, green, blue = split_into_rgb_channels(image)
    B = np.asmatrix(blue.astype(np.uint8))
    G = np.asmatrix(green.astype(np.uint8))
    R = np.asmatrix(red.astype(np.uint8))
    return B, G, R


def dna_encode(b, g, r):
    """Pixels -> bits -> DNA letters (2 bits per letter)."""
    b = np.unpackbits(b, axis=1)
    g = np.unpackbits(g, axis=1)
    r = np.unpackbits(r, axis=1)

    m, n = b.shape
    r_enc = np.empty((m, int(n / 2)), dtype="U1")
    g_enc = np.empty((m, int(n / 2)), dtype="U1")
    b_enc = np.empty((m, int(n / 2)), dtype="U1")

    for color, enc in zip((b, g, r), (b_enc, g_enc, r_enc)):
        idx = 0
        for j in range(0, m):
            for i in range(0, n, 2):
                enc[j, idx] = dna["{0}{1}".format(color[j, i], color[j, i + 1])]
                idx += 1
                if i == n - 2:
                    idx = 0
                    break

    return b_enc.astype(str), g_enc.astype(str), r_enc.astype(str)


def key_matrix_encode(key, b):
    """Build the key matrix from the hash and DNA-encode it."""
    b = np.unpackbits(b, axis=1)
    m, n = b.shape

    key_bin = bin(int(key, 16))[2:].zfill(256)

    Mk = np.zeros((m, n), dtype=np.uint8)
    x = 0
    for j in range(0, m):
        for i in range(0, n):
            Mk[j, i] = int(key_bin[x % 256])
            x += 1

    Mk_enc = np.empty((m, int(n / 2)), dtype="U1")
    idx = 0
    for j in range(0, m):
        for i in range(0, n, 2):
            if idx == (n / 2):
                idx = 0
            Mk_enc[j, idx] = dna["{0}{1}".format(Mk[j, i], Mk[j, i + 1])]
            idx += 1

    return Mk_enc.astype(str)


def xor_operation_encrypt(b, g, r, mk):
    """DNA XOR of the encoded channels with the encoded key matrix."""
    m, n = b.shape
    bx = np.empty((m, n), dtype="U1")
    gx = np.empty((m, n), dtype="U1")
    rx = np.empty((m, n), dtype="U1")

    b = b.astype(str)
    g = g.astype(str)
    r = r.astype(str)

    for i in range(0, m):
        for j in range(0, n):
            bx[i, j] = dna["{0}{1}".format(b[i, j], mk[i, j])]
            gx[i, j] = dna["{0}{1}".format(g[i, j], mk[i, j])]
            rx[i, j] = dna["{0}{1}".format(r[i, j], mk[i, j])]

    return bx.astype(str), gx.astype(str), rx.astype(str)


# ---------------------------------------------------------------------------
# Chaotic sequences
# ---------------------------------------------------------------------------

def gen_chaos_seq(m, n, x0, y0, z0):
    """Integrate the Lorenz system and return the x, y, z sequences."""
    N = m * n * 4
    t = np.linspace(0, TMAX, N)
    f = odeint(lorenz, (x0, y0, z0), t, args=(LORENZ_A, LORENZ_B, LORENZ_C))
    x, y, z = f.T
    return x[:N], y[:N], z[:N]


def plot(x, y, z):
    """Plot the Lorenz attractor (optional)."""
    import matplotlib.pyplot as plt

    N = len(x)
    plt.figure()
    ax = plt.axes(projection="3d")
    s = 100
    c = np.linspace(0, 1, N)
    for i in range(0, N - s, s):
        ax.plot(x[i:i + s + 1], y[i:i + s + 1], z[i:i + s + 1],
                color=(1 - c[i], c[i], 1), alpha=0.4)
    plt.show()


def sequence_indexing(x, y, z):
    """fx[i] = position of x[i] in the sorted order of x (same for y, z)."""
    n = len(x)
    fx = np.zeros(n, dtype=np.uint32)
    fy = np.zeros(n, dtype=np.uint32)
    fz = np.zeros(n, dtype=np.uint32)

    for seq, out in ((x, fx), (y, fy), (z, fz)):
        sorted_seq = sorted(seq)
        for k1 in range(0, n):
            out[k1] = bsearch(sorted_seq, seq[k1])

    return fx, fy, fz


# ---------------------------------------------------------------------------
# Scrambling
# ---------------------------------------------------------------------------

def scramble_encrypt(fx, fy, fz, b, r, g):
    """Scramble the blue, green and red DNA matrices using fz, fy and fx."""
    p, q = b.shape
    size = p * q

    bx = b.reshape(size).astype(str)
    gx = g.reshape(size).astype(str)
    rx = r.reshape(size).astype(str)

    bx_s = np.empty(size, dtype="U1")
    gx_s = np.empty(size, dtype="U1")
    rx_s = np.empty(size, dtype="U1")

    for i in range(size):
        bx_s[i] = bx[fz[i]]
    for i in range(size):
        gx_s[i] = gx[fy[i]]
    for i in range(size):
        rx_s[i] = rx[fx[i]]

    b_s = bx_s.astype(str).reshape(p, q)
    g_s = gx_s.astype(str).reshape(p, q)
    r_s = rx_s.astype(str).reshape(p, q)
    return b_s, g_s, r_s


def scramble_decrypt(fx, fy, fz, b, g, r):
    """Inverse of scramble_encrypt."""
    p, q = b.shape
    size = p * q

    bx = b.reshape(size).astype(str)
    gx = g.reshape(size).astype(str)
    rx = r.reshape(size).astype(str)

    bx_s = b.reshape(size).astype(str)
    gx_s = g.reshape(size).astype(str)
    rx_s = r.reshape(size).astype(str)

    for i in range(size):
        bx_s[fz[i]] = bx[i]
    for i in range(size):
        gx_s[fy[i]] = gx[i]
    for i in range(size):
        rx_s[fx[i]] = rx[i]

    return bx_s.reshape(p, q), gx_s.reshape(p, q), rx_s.reshape(p, q)


def dna_decode(b, g, r):
    """DNA letters -> bits -> pixel values."""
    m, n = b.shape
    r_dec = np.ndarray((m, int(n * 2)), dtype=np.uint8)
    g_dec = np.ndarray((m, int(n * 2)), dtype=np.uint8)
    b_dec = np.ndarray((m, int(n * 2)), dtype=np.uint8)

    for color, dec in zip((b, g, r), (b_dec, g_dec, r_dec)):
        for j in range(0, m):
            for i in range(0, n):
                dec[j, 2 * i] = dna["{0}".format(color[j, i])][0]
                dec[j, 2 * i + 1] = dna["{0}".format(color[j, i])][1]

    b_dec = np.packbits(b_dec, axis=-1)
    g_dec = np.packbits(g_dec, axis=-1)
    r_dec = np.packbits(r_dec, axis=-1)
    return b_dec, g_dec, r_dec


def xor_operation_decrypt(b, g, r, mk):
    """DNA XOR with the key matrix (XOR is its own inverse)."""
    m, n = b.shape
    bx = np.empty((m, n), dtype="U1")
    gx = np.empty((m, n), dtype="U1")
    rx = np.empty((m, n), dtype="U1")

    b = b.astype(str)
    g = g.astype(str)
    r = r.astype(str)

    for i in range(0, m):
        for j in range(0, n):
            bx[i, j] = dna["{0}{1}".format(b[i, j], mk[i, j])]
            gx[i, j] = dna["{0}{1}".format(g[i, j], mk[i, j])]
            rx[i, j] = dna["{0}{1}".format(r[i, j], mk[i, j])]

    return bx.astype(str), gx.astype(str), rx.astype(str)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _base_name(path):
    """File name without folders and extension: '/content/77.png' -> '77'."""
    return os.path.splitext(os.path.basename(path))[0]


def recover_image(b, g, r, iname):
    """Put the encrypted channels into the image and save it as PNG."""
    img = cv2.imread(iname)
    img[:, :, 2] = r
    img[:, :, 1] = g
    img[:, :, 0] = b
    out = "encrypt_{0}.png".format(_base_name(iname))
    cv2.imwrite(out, img)
    print("saved encrypted image as " + out)
    return img


def decrypt(image, fx, fy, fz, file_path, Mk, bt, gt, rt):
    """Decrypt the encrypted image and save the recovered one as PNG."""
    r, g, b = split_into_rgb_channels(image)

    p, q = rt.shape
    benc, genc, renc = dna_encode(b, g, r)
    bs, gs, rs = scramble_decrypt(fx, fy, fz, benc, genc, renc)
    bx, gx, rx = xor_operation_decrypt(bs, gs, rs, Mk)
    blue, green, red = dna_decode(bx, gx, rx)

    img = np.zeros((p, q, 3), dtype=np.uint8)
    img[:, :, 0] = blue
    img[:, :, 1] = green
    img[:, :, 2] = red

    out = "recovered_{0}.png".format(_base_name(file_path))
    cv2.imwrite(out, img)
    print("saved recovered image as " + out)
    return img


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(image_path=None, secret=None, show_plot=False):
    file_path = find_image(image_path) if image_path else image_input_from_user()
    secret = secret if secret is not None else key_input_from_user()

    # MD5 key + image size (m = width, n = height)
    key, m, n = securekey(file_path, secret)

    # Lorenz initial conditions from the key
    x0, y0, z0 = update_lorenz(key)

    # R, G, B matrices
    blue, green, red = decompose_matrix(file_path)

    # 1. DNA encoding of the channels
    blue_e, green_e, red_e = dna_encode(blue, green, red)

    # 2. DNA-encoded key matrix
    Mk_e = key_matrix_encode(key, blue)

    # 3. DNA XOR of the channels with the key matrix
    blue_x, green_x, red_x = xor_operation_encrypt(blue_e, green_e, red_e, Mk_e)

    # 4. Lorenz sequences and permutation indices
    x, y, z = gen_chaos_seq(m, n, x0, y0, z0)
    if show_plot:
        plot(x, y, z)
    fx, fy, fz = sequence_indexing(x, y, z)

    # 5. Scramble with the chaotic indices
    blue_s, green_s, red_s = scramble_encrypt(fx, fy, fz, blue_x, red_x, green_x)

    # 6. DNA decoding back to pixel values, then save the encrypted image
    b, g, r = dna_decode(blue_s, green_s, red_s)
    img = recover_image(b, g, r, file_path)

    # Decryption
    print("decrypting...")
    decrypt(img, fx, fy, fz, file_path, Mk_e, blue, green, red)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DNA + Lorenz + MD5 image encryption")
    parser.add_argument("--image", help="path to the input image")
    parser.add_argument("--key", help="secret key")
    parser.add_argument("--plot", action="store_true", help="show the Lorenz attractor")
    # parse_known_args: ignore extra arguments added by Jupyter/Colab (e.g. -f kernel.json)
    args, _ = parser.parse_known_args()
    main(args.image, args.key, args.plot)
