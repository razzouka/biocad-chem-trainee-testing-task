# Подготовка лигандов для AutoDock Vina

Решение тестового задания для стажировки в группу вычислительной химии BIOCAD.

## Описание

Скрипт `ligands_preparation.py` преобразует набор малых молекул лигандов из форматов SDF и SMILES в формат PDBQT для докинга в AutoDock Vina.
Для подготовки используются библиотеки RDKit (добавление водородов, генерация 3D-конформера, оптимизация UFF) и Meeko (создание PDBQT из RDKit-молекулы).

Поддерживаются SDF-файлы с 2D-координатами и без явных водородов; для каждого лиганда создаётся отдельный PDBQT-файл с именем лиганда.

## Установка окружения

Пример создания окружения (conda):

```bash
conda create -n ligand_prep -c conda-forge python=3.11
conda activate ligand_prep
conda install -c conda-forge rdkit numpy scipy gemmi
pip install meeko
```

Проверено на следующих версиях:
- Python 3.11.15
- RDKit 2026.03.1
- Meeko 0.7.1
- NumPy 2.4.3
- SciPy 1.17.1
- Gemmi 0.7.5

## Входные данные

В папке `raw_data/` должны быть файлы:

- `example.sdf` — SDF с несколькими лигандами  
- `example.smi` — SMILES-файл с лигандами (формат: SMILES и имя через пробел, первая строка — заголовок)

Примеры взяты из репозитория с тестовым заданием BIOCAD.

## Запуск

Примеры запуска:

```bash
# только SDF
python ligands_preparation.py --sdf raw_data/example.sdf

# только SMILES
python ligands_preparation.py --smi raw_data/example.smi

# SDF и SMILES
python ligands_preparation.py --sdf raw_data/example.sdf --smi raw_data/example.smi
```

### Параметры командной строки

- `--sdf <path>`  
  Путь к входному SDF‑файлу, содержащему структуры лигандов.

- `--smi <path>`  
  Путь к входному файлу SMILES (`.smi`), где каждая строка содержит одну SMILES‑запись молекулы.

Вы можете передать либо только `--sdf`, либо только `--smi`, либо оба параметра одновременно.  
Замените `path/to/file.sdf` и `path/to/file.smi` на реальные пути к вашим файлам (например, `raw_data/example.sdf`).

PDBQT-файлы сохраняются в:

- `sdf-pdbqt/` — для лигандов из SDF  
- `smi-pdbqt/` — для лигандов из SMILES  

Имена файлов соответствуют именам лигандов.

## Детали реализации

- RDKit: `Chem.SDMolSupplier` и `SmilesMolSupplier` для чтения SDF/SMILES, добавление водородов, 3D (ETKDG) и оптимизация UFF.  
- Meeko: `MoleculePreparation` и `PDBQTWriterLegacy` для подготовки лигандов и записи PDBQT, совместимых с AutoDock Vina и AutoDock-GPU. 
- Алифатические/макроциклические кольца с числом атомов больше 5 делаются гибкими (`rigid_macrocycles=False`, `min_ring_size=6`).
