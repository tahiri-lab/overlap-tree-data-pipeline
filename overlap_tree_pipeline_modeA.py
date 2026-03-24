#!/usr/bin/env python3
"""
Description:
    This script automates the creation of biologically meaningful partially overlapping phylogenetic trees
    with branch lengths using VertLife phylosubsets. Species selection supports three modes:
    user-provided list, uniform random sampling, and stratified random sampling by genus.
"""

import argparse
from collections import defaultdict
import math
import os
import random
import shutil
import sys
import time
import zipfile
from typing import Callable, Dict, List, Optional

import pandas as pd
import requests
from Bio import Phylo

# For Selenium
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service

# Argument parsing

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Construct partially overlapping phylogenetic tree datasets with branch lengths via VertLife."
    )
    parser.add_argument("species_group", choices=["amphibians", "birds", "mammals", "sharks", "squamates"])
    parser.add_argument("n", type=int, help="Number of species in the base set (ignored in user_list mode).")
    parser.add_argument("number_of_trees", type=int, help="Total number of trees in the final combined dataset.")
    parser.add_argument("email", type=str, help="Email for VertLife phylosubsets requests.")

    parser.add_argument(
        "--selection_mode",
        choices=["user_list", "random", "stratified"],
        default="stratified",
        help="Species selection mode used to construct the base set."
    )
    parser.add_argument(
        "--species_list_file",
        type=str,
        default=None,
        help="Path to a user-provided species list used when selection_mode=user_list."
    )
    parser.add_argument(
        "--stratify_by",
        choices=["genus"],
        default="genus",
        help="Taxonomic label used for stratified selection. Current implementation supports genus."
    )
    parser.add_argument(
        "--max_per_stratum",
        type=int,
        default=None,
        help="Optional cap on how many species can be drawn from any one stratum during stratified sampling."
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible selection and tree sampling."
    )
    return parser.parse_args()



# Overlap subset sizing helpers

def calculate_n_for_k(k: int, p_values: List[float]) -> int:
    n_total = k
    previous_new_species = 0
    for p in p_values:
        common_species = (2 * k * p) / (1 + p)
        rounded_common_species = round(common_species)
        new_species = k - rounded_common_species
        actual_new_species = new_species - previous_new_species
        n_total += math.ceil(actual_new_species)
        previous_new_species = new_species
    return int(n_total)


def find_k_from_n(n: int, max_k: int, p_values: List[float]) -> Optional[int]:
    for k in range(1, max_k + 1):
        if calculate_n_for_k(k, p_values) >= n:
            return k
    return None


# Species selection

def genus_label(species_name: str) -> str:
    parts = species_name.split()
    return parts[0] if parts else species_name


def load_species_list_file(path: str) -> List[str]:
    if not path or not os.path.exists(path):
        raise FileNotFoundError("species_list_file was not provided or the file does not exist.")

    _, ext = os.path.splitext(path.lower())

    if ext in [".txt", ".list"]:
        with open(path, "r", encoding="utf-8") as f:
            return [ln.strip() for ln in f if ln.strip()]

    if ext == ".csv":
        df = pd.read_csv(path)
        if df.shape[1] == 0:
            return []
        col0 = df.columns[0]
        return df[col0].dropna().astype(str).tolist()

    with open(path, "r", encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]


def stratified_random_sample(
    items: List[str],
    n: int,
    label_fn: Callable[[str], str],
    max_per_stratum: Optional[int] = None
) -> List[str]:
    if n > len(items):
        raise ValueError(f"Requested n={n} exceeds the available species count {len(items)}.")

    strata: Dict[str, List[str]] = defaultdict(list)
    for sp in items:
        strata[label_fn(sp)].append(sp)

    labels = list(strata.keys())
    for lbl in labels:
        random.shuffle(strata[lbl])
    random.shuffle(labels)

    selected: List[str] = []
    picked_per: Dict[str, int] = defaultdict(int)

    active = labels[:]  # labels that still have available species
    idx = 0
    while len(selected) < n and active:
        lbl = active[idx % len(active)]

        if not strata[lbl]:
            active.remove(lbl)
            continue

        if max_per_stratum is not None and picked_per[lbl] >= max_per_stratum:
            idx += 1
            # Stop if every remaining stratum is blocked by max_per_stratum or empty
            blocked = True
            for x in active:
                if strata[x] and (max_per_stratum is None or picked_per[x] < max_per_stratum):
                    blocked = False
                    break
            if blocked:
                break
            continue

        selected.append(strata[lbl].pop())
        picked_per[lbl] += 1
        idx += 1

    if len(selected) < n:
        raise ValueError(
            "Stratified sampling did not reach the requested n under the current max_per_stratum constraint."
        )
    return selected


def select_base_species(
    full_species_list: List[str],
    n: int,
    selection_mode: str,
    species_list_file: Optional[str] = None,
    stratify_by: str = "genus",
    max_per_stratum: Optional[int] = None
) -> List[str]:
    if selection_mode == "user_list":
        user_species = load_species_list_file(species_list_file)
        user_species = [s.strip() for s in user_species if isinstance(s, str) and s.strip()]
        user_species = list(dict.fromkeys(user_species))  # remove duplicates, preserve order

        unknown = [s for s in user_species if s not in full_species_list]
        if unknown:
            preview = ", ".join(unknown[:10])
            raise ValueError(f"User list contains taxa not present in the group list. Examples: {preview}")

        return user_species

    if selection_mode == "random":
        if n > len(full_species_list):
            raise ValueError(f"Requested n={n} exceeds the available species count {len(full_species_list)}.")
        return random.sample(full_species_list, n)

    if selection_mode == "stratified":
        if stratify_by != "genus":
            raise ValueError(f"Unsupported stratify_by option: {stratify_by}")
        return stratified_random_sample(full_species_list, n, genus_label, max_per_stratum=max_per_stratum)

    raise ValueError(f"Unknown selection_mode: {selection_mode}")


# Overlapping subset construction

def create_overlapping_subsets(
    df: pd.DataFrame,
    p_values: List[float]
) -> None:
    for group in df.columns:
        species_list = df[group].dropna().tolist()
        n_local = len(species_list)

        k = find_k_from_n(n_local, max_k=1000, p_values=p_values)
        if k is None:
            print(f"Could not find a valid subset size for group {group}")
            continue

        subsets: Dict[str, List[Optional[str]]] = {f"Subset {i}": [] for i in range(1, 11)}
        start_index = 0

        subsets["Subset 1"] = species_list[start_index:start_index + k]

        for i, p in enumerate(p_values, start=2):
            common_species = round((2 * k * p) / (1 + p))
            new_species = k - common_species

            subset_start = start_index + int(new_species)
            subset_end = subset_start + k

            subset_species = species_list[subset_start:subset_end]
            if len(subset_species) < k:
                subset_species += species_list[:k - len(subset_species)]

            subsets[f"Subset {i}"] = subset_species

        max_len = len(species_list)
        for key in subsets:
            subsets[key] += [None] * (max_len - len(subsets[key]))

        output_df = pd.DataFrame(subsets)
        group_slug = str(group).strip().lower()
        output_file = f"{group_slug}_overlapping_subsets.csv"
        output_df.to_csv(output_file, index=False)
        print(f"Saved file for {group}: {output_file}")


# VertLife interaction

GROUP_MAPPING = {
    "amphibians": "amphibiantree",
    "birds": "birdtree",
    "mammals": "mammaltree",
    "sharks": "sharktree",
    "squamates": "squamatetree"
}

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
import os

def setup_driver() -> webdriver.Chrome:
    chrome_options = Options()
    chrome_options.add_argument("--headless=new")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--window-size=1920x1080")

    # Linux flags
    if os.name != "nt":
        chrome_options.add_argument("--no-sandbox")
        chrome_options.add_argument("--disable-gpu")
        chrome_options.add_argument("--disable-setuid-sandbox")
        chrome_options.add_argument("--remote-debugging-port=9222")

    # Optional explicit browser path
    chrome_bin = os.environ.get("CHROME_BIN")
    if chrome_bin and os.path.exists(chrome_bin):
        chrome_options.binary_location = chrome_bin

    # Optional explicit driver path
    chromedriver_bin = os.environ.get("CHROMEDRIVER")
    if chromedriver_bin and os.path.exists(chromedriver_bin):
        service = Service(executable_path=chromedriver_bin)
        return webdriver.Chrome(service=service, options=chrome_options)

    # Best default: let Selenium Manager resolve/download what it needs
    return webdriver.Chrome(options=chrome_options)

def clean_folder(folder_name: str) -> None:
    if os.path.exists(folder_name):
        print(f"Cleaning folder: {folder_name}")
        shutil.rmtree(folder_name)
    os.makedirs(folder_name, exist_ok=True)


def submit_tree_request(driver: webdriver.Chrome, species_list: List[str], email: str, group: str) -> Optional[str]:
    driver.get("https://vertlife.org/phylosubsets/")

    try:
        if group.lower() != "amphibians":
            print(f"Selecting species group: {group}")
            group_value = GROUP_MAPPING.get(group.lower())
            if not group_value:
                print(f"Group {group} not found in the mapping.")
                return None

            group_radio_button = WebDriverWait(driver, 30).until(
                EC.element_to_be_clickable((By.XPATH, f"//input[@type='radio' and @value='{group_value}']"))
            )
            group_radio_button.click()
        else:
            print("Skipping species group selection since amphibians is the default.")

        species_textarea = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.ID, "selected"))
        )
        species_textarea.clear()
        species_textarea.send_keys("\n".join(species_list))

        print(f"Email being used: {email}")
        email_input = driver.find_element(By.ID, "email")
        email_input.clear()
        email_input.send_keys(email)

        print("Submitting tree request...")
        get_trees_button = driver.find_element(By.ID, "btnGetTrees")
        get_trees_button.click()

        print("Waiting for job ID...")
        job_id_element = WebDriverWait(driver, 60).until(
            EC.presence_of_element_located((By.XPATH, "//div[@id='status']//strong"))
        )
        job_id = job_id_element.text.strip()
        print(f"Job ID: {job_id}")
        return job_id

    except Exception as e:
        print(f"Error during tree request: {e}")
        return None


def download_zipfile_using_job_id(job_id: str, folder_name: str, max_retries: int = 20) -> bool:
    download_url = f"https://data.vertlife.org/pruned_treesets/{job_id}/{job_id}.zip"
    print(f"Downloading from: {download_url}")

    retries = 0
    while retries < max_retries:
        response = requests.get(download_url, timeout=120)

        if response.status_code == 200:
            zip_filename = os.path.join(folder_name, f"{job_id}.zip")
            with open(zip_filename, "wb") as f:
                f.write(response.content)

            try:
                with zipfile.ZipFile(zip_filename, "r") as zip_ref:
                    print(f"Extracting files from {zip_filename}...")
                    for file in zip_ref.namelist():
                        file_extension = os.path.splitext(file)[1]
                        new_file_name = f"{job_id}{file_extension}"
                        target_path = os.path.join(folder_name, new_file_name)

                        zip_ref.extract(file, folder_name)
                        os.rename(os.path.join(folder_name, file), target_path)
                        print(f"Extracted: {target_path}")
            except zipfile.BadZipFile:
                print(f"Error: {zip_filename} is not a valid zip file.")
                return False

            return True

        retries += 1
        print(
            f"Failed to download from {download_url}. Status code: {response.status_code}. "
            f"Retrying ({retries}/{max_retries})"
        )
        time.sleep(30)

    print(f"Failed to download file after {max_retries} attempts.")
    return False


def automate_process(input_file: str, email: str) -> None:
    group_name = os.path.basename(input_file).split("_")[0].strip().lower()
    print(f"Group Name: {group_name}")

    folder_name = f"{group_name}_nexus"
    clean_folder(folder_name)

    df = pd.read_csv(input_file)
    job_ids: Dict[str, str] = {}

    for subset_name in df.columns:
        success = False
        retries = 0
        while not success and retries < 15:
            try:
                print(f"Requesting trees for {subset_name}")
                subset_species = df[subset_name].dropna().tolist()

                driver = setup_driver()
                job_id = submit_tree_request(driver, subset_species, email, group_name)
                driver.quit()

                if job_id:
                    job_ids[subset_name] = job_id
                    success = True
                else:
                    print(f"Retrying {subset_name} due to missing job ID.")
                    retries += 1
                    time.sleep(20)

            except Exception as e:
                print(f"Exception during processing subset {subset_name}: {e}")
                retries += 1
                time.sleep(20)

        time.sleep(10)

    print("Waiting for 60 seconds before starting downloads to ensure files are ready...")
    time.sleep(60)

    for subset_name, job_id in job_ids.items():
        print(f"Processing download for {subset_name}")
        ok = download_zipfile_using_job_id(job_id, folder_name)
        if not ok:
            print(f"Failed to download trees for {subset_name}")


# Nexus to Newick conversion and combining

def convert_nexus_to_newick(nexus_file: str, t: int) -> List[str]:
    try:
        trees = list(Phylo.parse(nexus_file, "nexus"))
    except Exception as e:
        print(f"Error parsing {nexus_file}: {e}")
        return []

    if len(trees) < t:
        print(f"Warning: {nexus_file} has fewer than {t} trees. Selecting all available trees.")
        t = len(trees)

    selected_trees = random.sample(trees, t)

    selected_newick_trees: List[str] = []
    for tree in selected_trees:
        newick_str = tree.format("newick").strip()
        if newick_str.endswith(":0.00000;"):
            newick_str = newick_str[:-9] + ";"
        selected_newick_trees.append(newick_str)

    return selected_newick_trees


def process_nexus_files(input_dir: str, output_file: str, t: int) -> None:
    all_selected_trees: List[str] = []

    for filename in os.listdir(input_dir):
        if filename.endswith(".nex"):
            nexus_file = os.path.join(input_dir, filename)
            selected_trees = convert_nexus_to_newick(nexus_file, t)
            all_selected_trees.extend(selected_trees)

    with open(output_file, "w", encoding="utf-8") as out_file:
        for tree in all_selected_trees:
            if tree:
                out_file.write(tree + "\n")

    print(f"Saved {len(all_selected_trees)} Newick trees to {output_file}")


# Main

def main() -> None:
    args = parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    species_group = args.species_group.lower()
    n = int(args.n)
    number_of_trees = int(args.number_of_trees)
    email = args.email

    if number_of_trees < 10 or number_of_trees > 1000 or number_of_trees % 10 != 0:
        print("Error: The number of trees must be between 10 and 1000 and divisible by 10.")
        sys.exit(1)

    t = number_of_trees // 10  # number of trees to select per subset

    if not os.path.exists("all_species_lists.csv"):
        print("Error: 'all_species_lists.csv' file not found. Please ensure the file is present in the current directory.")
        sys.exit(1)

    data = pd.read_csv("all_species_lists.csv")

    species_dict = {
        "amphibians": data["Amphibians"].dropna().tolist(),
        "birds": data["Birds"].dropna().tolist(),
        "mammals": data["Mammals"].dropna().tolist(),
        "sharks": data["Sharks"].dropna().tolist(),
        "squamates": data["Squamates"].dropna().tolist(),
    }

    species_list = species_dict[species_group]

    try:
        final_species = select_base_species(
            species_list,
            n,
            args.selection_mode,
            species_list_file=args.species_list_file,
            stratify_by=args.stratify_by,
            max_per_stratum=args.max_per_stratum
        )
        if args.selection_mode == "user_list":
            n = len(final_species)
    except Exception as e:
        print(f"Error during species selection: {e}")
        sys.exit(1)

    result_df = pd.DataFrame({species_group: final_species})
    result_df.to_csv("selected_species.csv", index=False)

    p_values = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    df = pd.read_csv("selected_species.csv")
    create_overlapping_subsets(df, p_values=p_values)

    input_file = f"{species_group}_overlapping_subsets.csv"
    automate_process(input_file, email)

    input_dir = f"./{species_group}_nexus/"
    output_file = f"overlapping_dataset_{species_group}.txt"
    process_nexus_files(input_dir, output_file, t)

    citations = {
        "amphibians": "Amphibians: Jetz, W., & Pyron, R. A. (2018). The interplay of past diversification and evolutionary isolation with present imperilment across the amphibian tree of life. Nature ecology & evolution, 2(5), 850-858.",
        "birds": "Birds: Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The global diversity of birds in space and time. Nature, 491(7424), 444-448.",
        "mammals": "Mammals: Upham, N. S., Esselstyn, J. A., & Jetz, W. (2019). Inferring the mammal tree: species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS biology, 17(12), e3000494.",
        "sharks": "Sharks: Stein, R. W., Mull, C. G., Kuhn, T. S., Aschliman, N. C., Davidson, L. N., Joy, J. B., ... & Mooers, A. O. (2018). Global priorities for conserving the evolutionary history of sharks, rays and chimaeras. Nature ecology & evolution, 2(2), 288-298.",
        "squamates": "Squamates: Tonini, J. F. R., Beard, K. H., Ferreira, R. B., Jetz, W., & Pyron, R. A. (2016). Fully-sampled phylogenies of squamates reveal evolutionary patterns in threat status. Biological Conservation, 204, 23-31.",
    }

    print(f"Dataset saved to {output_file}")
    print("Please cite the following source for the data used:")
    print(citations[species_group])
    print("Website with comprehensive data: https://vertlife.org/data/")


if __name__ == "__main__":
    main()
