from bs4 import BeautifulSoup
import requests
import fnmatch
import zipfile
import io
import os


def get_page(url):
    page = requests.get(url)
    parsed = BeautifulSoup(page.content, 'html.parser')

    return parsed


def save_version(version, url, save_to):
    version = version.replace('.', '_')

    if '-2_' in version:
        "\tIgnoring - no tab2 file as version 2"
        return

    print(f"getting {version}")

    version = save_to

    try:
        os.listdir(version)
    except FileNotFoundError:
        os.mkdir(version)

    if os.listdir(version):
        return

    directory = get_page(url)
    downloads = directory.find_all('a', class_="filePopup")

    download = fnmatch.filter(map(lambda x: x['href'], downloads), '*ORGANISM*.tab2.zip')

    if download:
        download = download[0]
    else:
        print("\tNo file in tab2 format found")
        return

    print("\tDownloading file")
    r = requests.get(download, stream=True)

    if r.ok:
        print("\tUnzipping file")
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(version)
    else:
        print("\tDownload failed")
        return


if __name__ == "__main__":
    archive_url = "https://downloads.thebiogrid.org/BioGRID/Release-Archive/"
    archive_directory = "biogrid_archive/"

    page_html = get_page(archive_url)

    links = page_html.find_all('a', title="View Contents of this Directory")

    link_dict = {link.text: link['href'] for link in links}

    for key, value in sorted(link_dict.items()):
        try:
            save_version(key, value, save_to=archive_directory)
        except zipfile.BadZipfile:
            print("\tError unzipping file - try again")
            save_version(key, value, save_to=archive_directory)

    print("Done!")
