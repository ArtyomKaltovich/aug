import urllib.request
import json
import os
import pandas as pd

from aug.data.fasta import parse_fasta_string


def request_protein_sequence(*ids, source: str = "uniprot"):
    assert source == "uniprot", "currently other sources isn't supported, you can use uniprot only"
    for id in ids:
        content = urllib.request.urlopen(f"https://www.uniprot.org/uniprot/{id}.fasta").read().decode()
        yield from parse_fasta_string(content)


def request_and_parse_json(url, json_callback=None, **kwargs):
    """ Request a web-page with json content and parse it.
    It can be usable for accessing web api.
    :param url: an url to request
    :param kwargs: parameter of query, they will be transform into GET request params (e.g. name1=value1&name2=value2)
    :return: parsed JSON as a dict object
    """
    params = [str(key) + "=" + str(value) for key, value in kwargs.items()]
    params = "&".join(params)
    content = urllib.request.urlopen(f"{url}?{params}").read().decode()
    if json_callback:
        content = json_callback(content)
    content = json.loads(content)
    return content


def download_from_mg_rast(folder_to_save_: str, just_description_=False, files_to_download_=None,
                          description_callback_=None, rewrite_existed_=False, path_="", **kwargs):
    """

    :param folder_to_save_:
    :param just_description_:
    :param files_to_download_:
    :param description_callback_:
    :param rewrite_existed_: rewrite file if file or directory with such name is already presented
    :param path_: write in the middle of path (can be used for folder creating)
    :param kwargs:
    :return:
    """
    def _save_json_to_file_callback(content):
        with open(folder_to_save_ + "metagenomes_descr.json", "w") as f:
            f.write(content)
        return content

    assert files_to_download_, "Please, specify files types which should be downloaded"

    folder_to_save_ = folder_to_save_ if folder_to_save_.endswith(os.sep) else folder_to_save_ + os.sep
    metagenomes_descr = request_and_parse_json(url="http://api.mg-rast.org/search",
                                               json_callback=_save_json_to_file_callback,
                                               **kwargs)
    metagenomes_descr = metagenomes_descr["data"]
    _save_mg_rast_descr(metagenomes_descr, folder_to_save_ + "metagenomes_descr.csv", description_callback_,
                        rewrite_existed_)

    if not just_description_:
        _download_mg_rast_files(files_to_download_, folder_to_save_, metagenomes_descr, path_, rewrite_existed_)


def _download_mg_rast_files(files_to_download_, folder_to_save_, metagenomes_descr, path_, rewrite_existed_):
    with open(folder_to_save_ + "metagenomes_details_urls.csv", "w") as csv:
        csv.write("id, path, details_url\n")
        for i, record in enumerate(metagenomes_descr, 1):
            try:
                _process_mg_rast_metagenome_record(csv, files_to_download_, folder_to_save_, record,
                                                   rewrite_existed_, path_)
                print(f"{i}/{len(metagenomes_descr)}")
            except urllib.error.HTTPError as e:
                print(e)


def _process_mg_rast_metagenome_record(csv, files_to_download_, folder_to_save_, record_, rewrite_existed_, path_):
    id = record_["metagenome_id"]
    record_info = request_and_parse_json(url=f"http://api.mg-rast.org/download/{id}")
    folder_to_save_ = _update_folder_to_save(folder_to_save_, path_, record_)
    _download_mg_rast_files_by_postfixes(folder_to_save_, record_info, files_to_download_, rewrite_existed_)
    csv.write(f"{id}, https://www.mg-rast.org/mgmain.html?mgpage=overview&metagenome={id}\n")


def _update_folder_to_save(folder_to_save_, path_, record):
    if path_:
        for p in path_:
            folder_to_save_ += record[p] + os.sep
            if not os.path.exists(folder_to_save_):
                os.mkdir(folder_to_save_)
    return folder_to_save_


def _save_mg_rast_descr(metagenomes_descr, path_to_save, callback, rewrite_existed_):
    data = pd.DataFrame(metagenomes_descr)
    if callback:
        callback(data)
    mode = "w" if rewrite_existed_ else "a"
    data.to_csv(path_to_save, mode=mode)


def _download_mg_rast_files_by_postfixes(folder_to_save_, record_info_, files_to_download_, rewrite_existed_):
    for file in record_info_["data"]:
        for type_ in files_to_download_:
            file_name = file["file_name"]
            if file["file_name"].endswith(type_):
                url_to_file = file["url"]
                file_path = folder_to_save_ + file_name
                if not os.path.exists(file_path) or rewrite_existed_:
                    urllib.request.urlretrieve(url_to_file, file_path)
                break
