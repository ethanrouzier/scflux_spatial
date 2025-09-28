# Tests pour scflux_spatial

Ce répertoire contient tous les tests unitaires et d'intégration pour le package `scflux_spatial`.

## Structure des tests

```
tests/
├── __init__.py              # Initialisation du module de tests
├── conftest.py              # Configuration pytest et fixtures
├── test_gpr.py              # Tests pour le parsing GPR
├── test_fba.py              # Tests pour l'analyse de flux (FBA)
├── test_rd.py               # Tests pour les équations réaction-diffusion
├── test_performance.py      # Tests de performance et scalabilité
├── test_data/               # Données de test (si nécessaire)
└── README.md                # Cette documentation
```

## Types de tests

### 1. Tests unitaires (`test_*.py`)

- **`test_gpr.py`** : Tests pour le parsing et l'évaluation des règles GPR
  - Opérateurs AND/OR
  - Expressions imbriquées complexes
  - Gestion des erreurs et cas limites
  - Opérateurs personnalisés

- **`test_fba.py`** : Tests pour l'analyse de flux métaboliques
  - Modèle toy avec solution analytique connue
  - pFBA (minimisation de la somme des flux absolus)
  - Intégration de l'expression génique (E-Flux, iMAT)
  - Validation des contraintes et bornes

- **`test_rd.py`** : Tests pour les équations réaction-diffusion
  - Diffusion pure avec solutions analytiques connues
  - Conditions aux limites (Dirichlet, Neumann)
  - Convergence vers l'état stationnaire
  - Conservation de la masse

- **`test_performance.py`** : Tests de performance et scalabilité
  - Temps d'exécution avec grandes bases de données
  - Utilisation mémoire
  - Évaluation concurrente
  - Scalabilité avec la taille des données

### 2. Configuration des tests (`conftest.py`)

Contient les fixtures pytest communes :
- `mock_gene_expression` : Données d'expression génique simulées
- `mock_gpr_rules` : Règles GPR de test
- `mock_flux_bounds` : Bornes de flux simulées
- `mock_toy_model` : Modèle métabolique toy
- `mock_spatial_coordinates` : Coordonnées spatiales simulées
- `mock_concentration_field` : Champ de concentration simulé
- `mock_rd_parameters` : Paramètres RD simulés

### 3. Marquage des tests

Les tests sont marqués avec des étiquettes pour faciliter l'exécution sélective :

- `@pytest.mark.unit` : Tests unitaires (rapides, isolés)
- `@pytest.mark.integration` : Tests d'intégration (plus lents)
- `@pytest.mark.slow` : Tests lents (plusieurs secondes)
- `@pytest.mark.performance` : Tests de performance

## Exécution des tests

### 1. Tous les tests

```bash
# Avec pytest
python -m pytest tests/

# Avec le script de test
python run_tests.py

# Tests unitaires uniquement
python run_tests.py --type unit

# Tests rapides uniquement (exclure les tests lents)
python run_tests.py --type fast
```

### 2. Tests spécifiques

```bash
# Test spécifique
python -m pytest tests/test_gpr.py

# Test avec marqueur spécifique
python -m pytest -m "unit"

# Test avec couverture
python -m pytest --cov=scflux_spatial --cov-report=html
```

### 3. Avec unittest

```bash
# Test spécifique avec unittest
python -m unittest tests.test_gpr -v

# Tous les tests avec unittest
python -m unittest discover tests/ -v
```

### 4. Tests de performance

```bash
# Tests de performance uniquement
python -m pytest -m "performance" tests/

# Tests lents uniquement
python run_tests.py --type slow
```

## Configuration pytest

Le fichier `pytest.ini` configure :
- Découverte automatique des tests
- Marqueurs personnalisés
- Options de sortie
- Filtrage des avertissements
- Configuration de couverture

## Fixtures et données de test

### Fixtures principales

```python
@pytest.fixture
def mock_gene_expression():
    """Données d'expression génique simulées."""
    return {
        'GENE_A': 5.0,
        'GENE_B': 3.0,
        'GENE_C': 2.0,
        # ...
    }

@pytest.fixture
def mock_toy_model():
    """Modèle métabolique toy pour les tests."""
    # Modèle simple avec solution analytique connue
    pass
```

### Utilisation des fixtures

```python
def test_gpr_evaluation(mock_gene_expression, mock_gpr_rules):
    """Test d'évaluation GPR avec fixtures."""
    gpr_parser = GPRParser()
    
    result = gpr_parser._evaluate_gpr_rule_with_operators(
        mock_gpr_rules['simple_and'],
        mock_gene_expression
    )
    
    assert result == 3.0  # min(5.0, 3.0)
```

## Tests de validation

### 1. Tests GPR

- **Opérateurs logiques** : AND (min), OR (max)
- **Expressions imbriquées** : `(A and B) or (C and D)`
- **Cas limites** : valeurs nulles, négatives, très grandes
- **Gestion d'erreurs** : gènes manquants, expressions malformées
- **Opérateurs personnalisés** : fonctions AND/OR personnalisées

### 2. Tests FBA

- **Modèle toy** : Chaîne linéaire A → B → C → D
- **Solution analytique** : Tous les flux égaux, limités par le goulot
- **pFBA** : Minimisation de Σ|v| avec objectif maintenu
- **Intégration expression** : E-Flux, iMAT, méthodes linéaires
- **Validation physique** : Conservation de masse, contraintes

### 3. Tests RD

- **Diffusion pure** : Solution uniforme avec conditions Dirichlet
- **Conservation de masse** : Conditions Neumann
- **Convergence** : Critère L2 pour l'état stationnaire
- **Solutions analytiques** : Disque, géométries simples
- **Paramètres physiques** : Unités cohérentes, valeurs réalistes

## Métriques de performance

### 1. Temps d'exécution

- **GPR** : < 1 seconde pour 1000 gènes
- **FBA** : < 1 seconde pour modèle toy
- **RD** : < 10 secondes pour grille 200×200

### 2. Utilisation mémoire

- **Augmentation** : < 500 MB pour datasets de 10000 éléments
- **Efficacité** : Opérations vectorisées NumPy
- **Scalabilité** : Linéaire avec la taille des données

### 3. Couverture de code

- **Objectif** : > 90% de couverture
- **Rapport** : HTML généré avec `--cov-report=html`
- **Exclusions** : Tests, fichiers de configuration

## Débogage des tests

### 1. Tests qui échouent

```bash
# Mode verbose pour plus de détails
python -m pytest tests/test_gpr.py -v

# Arrêter au premier échec
python -m pytest tests/ -x

# Afficher les valeurs des variables
python -m pytest tests/ --tb=long
```

### 2. Tests lents

```bash
# Identifier les tests lents
python -m pytest tests/ --durations=10

# Exclure les tests lents
python -m pytest tests/ -m "not slow"
```

### 3. Problèmes de mémoire

```bash
# Tests de performance avec monitoring mémoire
python -m pytest tests/test_performance.py -v
```

## Contribution aux tests

### 1. Ajouter de nouveaux tests

1. Créer un nouveau fichier `test_*.py`
2. Importer les fixtures nécessaires
3. Utiliser les marqueurs appropriés
4. Ajouter la documentation

### 2. Bonnes pratiques

- **Nommage** : `test_*` pour les fonctions, `Test*` pour les classes
- **Isolation** : Chaque test doit être indépendant
- **Fixtures** : Réutiliser les fixtures communes
- **Assertions** : Messages d'erreur informatifs
- **Performance** : Marquer les tests lents

### 3. Exemple de test

```python
import pytest
from scflux_spatial.gem.gpr import GPRParser

class TestNewFeature(unittest.TestCase):
    """Tests pour une nouvelle fonctionnalité."""
    
    def setUp(self):
        """Configuration des tests."""
        self.gpr_parser = GPRParser()
    
    @pytest.mark.unit
    def test_new_functionality(self, mock_gene_expression):
        """Test de la nouvelle fonctionnalité."""
        # Arrange
        expected_result = 42
        
        # Act
        actual_result = self.gpr_parser.new_method(mock_gene_expression)
        
        # Assert
        self.assertEqual(actual_result, expected_result)
```

## Ressources

- [Documentation pytest](https://docs.pytest.org/)
- [Fixtures pytest](https://docs.pytest.org/en/stable/fixture.html)
- [Marqueurs pytest](https://docs.pytest.org/en/stable/mark.html)
- [Couverture de code](https://coverage.readthedocs.io/)
