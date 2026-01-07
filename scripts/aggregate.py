#!/usr/bin/env python3
"""
Aggregation and Leaderboard Generator
======================================
Collects all evaluation results and generates summary reports.

Outputs:
- CSV leaderboard with all scores
- HTML interactive report
- JSON summary with statistics
"""

import argparse
import json
import logging
from collections import defaultdict
from datetime import datetime
from pathlib import Path


try:
    import pandas as pd

    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def load_eval_results(eval_dir: Path) -> list[dict]:
    """Load all evaluation result JSON files and extract multiple metrics.

    Expected directory structure: {eval_dir}/{featurizer}/{predictor}/{target_id}/{scorer}.json
    """
    results = []

    # Dictionary to map file scorer names to internal metric keys
    # If the JSON file contains the metric directly, we use it.
    # Otherwise we might need to map specific file contents.

    for json_file in eval_dir.rglob("*.json"):
        try:
            with open(json_file) as f:
                data = json.load(f)

            # Extract info from path: {eval_dir}/{featurizer}/{predictor}/{target_id}/{scorer}.json
            parts = json_file.relative_to(eval_dir).parts

            # Check for the expected 4-level structure
            if len(parts) == 4:
                featurizer = parts[0]
                predictor = parts[1]
                target_id = parts[2]
                scorer_file = parts[3]  # e.g. "ost.json" or "dockq.json"
                scorer_name = json_file.stem  # e.g. "ost" or "dockq"

                method = f"{predictor}_{featurizer}"
            else:
                # Fallback for old/flat structure (though likely not used given Snakefile)
                logger.warning(
                    f"Unexpected directory structure for {json_file}, skipping parsing path metadata."
                )
                continue

            # If the file is a specific scorer file (e.g. ost.json), it might contain multiple metrics
            # We treat the file stem as the primary 'scorer' category, but we can extract multiple
            # specific metrics from it.

            # Common metrics we want to extract if present
            # For header generation later, we will use "scorer_name" from here.
            # But the 'scorer' column in the CSV is actually the specific metric (e.g. 'lddt', 'tm_score')

            # Let's look at what's in the data.
            # If data is a dict, iterate keys.

            if isinstance(data, dict):
                # Helper to add a score result
                def add_result(metric_name, metric_value):
                    # Handle lists (mean)
                    if isinstance(metric_value, list):
                        try:
                            metric_value = (
                                sum(metric_value) / len(metric_value)
                                if metric_value
                                else 0.0
                            )
                        except TypeError:
                            return  # Skip non-numeric

                    results.append(
                        {
                            "target_id": target_id,
                            "method": method,
                            "predictor": predictor,
                            "featurizer": featurizer,
                            "scorer_name": metric_name,
                            "score": float(metric_value),
                            "_source_file": str(json_file),
                        }
                    )

                # 1. Top level keys
                for key, value in data.items():
                    if isinstance(value, (int, float, list)) and key not in [
                        "score_details"
                    ]:
                        # Avoid duplicating if 'score' is just a summary
                        # But user wants everything so let's keep it unless it's redundant?
                        add_result(key, value)

                # 2. Nested score_details
                if "score_details" in data and isinstance(data["score_details"], dict):
                    for key, value in data["score_details"].items():
                        if isinstance(value, (int, float, list)):
                            add_result(key, value)

        except json.JSONDecodeError as e:
            logger.warning(f"Failed to parse {json_file}: {e}")
        except Exception as e:
            logger.warning(f"Error processing {json_file}: {e}")

    logger.info(f"Loaded {len(results)} evaluation results from {eval_dir}")
    return results


def build_leaderboard(results: list[dict], config: dict) -> dict:
    """Build leaderboard structure from results."""
    # Group by method (predictor_featurizer)
    method_scores = defaultdict(lambda: defaultdict(list))

    for result in results:
        method = result.get("method", "unknown")
        scorer = result.get("scorer_name", "unknown")
        score = result.get("score", 0.0)
        target = result.get("target_id", "unknown")

        method_scores[method][scorer].append(
            {
                "target_id": target,
                "score": score,
            }
        )

    # Compute aggregated statistics
    leaderboard = {
        "generated_at": datetime.now().isoformat(),
        "config": {
            "enabled_models": config.get("enabled_models", []),
            "enabled_scorers": config.get("enabled_scorers", []),
        },
        "methods": {},
    }

    for method, scorer_results in method_scores.items():
        method_stats = {
            "total_targets": 0,
            "scorers": {},
        }

        for scorer, scores in scorer_results.items():
            score_values = [s["score"] for s in scores]

            method_stats["scorers"][scorer] = {
                "mean": sum(score_values) / len(score_values) if score_values else 0,
                "median": sorted(score_values)[len(score_values) // 2]
                if score_values
                else 0,
                "min": min(score_values) if score_values else 0,
                "max": max(score_values) if score_values else 0,
                "count": len(score_values),
                "success_rate": sum(1 for s in score_values if s > 0.5)
                / len(score_values)
                if score_values
                else 0,
            }
            method_stats["total_targets"] = max(
                method_stats["total_targets"], len(scores)
            )

        leaderboard["methods"][method] = method_stats

    return leaderboard


def generate_csv(results: list[dict], output_path: Path):
    """Generate CSV leaderboard."""
    if HAS_PANDAS:
        df = pd.DataFrame(results)

        # Pivot to wide format
        if not df.empty:
            pivot = df.pivot_table(
                index=["target_id", "predictor", "featurizer"],
                columns="scorer_name",
                values="score",
                aggfunc="first",
            ).reset_index()
            pivot.to_csv(output_path, index=False)
        else:
            # Empty file
            with open(output_path, "w") as f:
                f.write("target_id,predictor,featurizer\n")
    else:
        # Manual CSV generation
        with open(output_path, "w") as f:
            if results:
                # Get all unique scorers
                scorers = sorted(set(r.get("scorer_name", "") for r in results))

                # Header
                f.write("target_id,predictor,featurizer," + ",".join(scorers) + "\n")

                # Group by target/predictor/featurizer
                grouped = defaultdict(dict)
                for r in results:
                    key = (
                        r.get("target_id", ""),
                        r.get("predictor", ""),
                        r.get("featurizer", ""),
                    )
                    grouped[key][r.get("scorer_name", "")] = r.get("score", 0)

                # Rows
                for (target, predictor, featurizer), scores in sorted(grouped.items()):
                    row = [target, predictor, featurizer] + [
                        str(scores.get(s, "")) for s in scorers
                    ]
                    f.write(",".join(row) + "\n")

    logger.info(f"Generated CSV: {output_path}")


def generate_html(leaderboard: dict, output_path: Path):
    """Generate HTML report with interactive tables."""
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PSP Benchmark Leaderboard</title>
    <style>
        :root {{
            --bg-primary: #0d1117;
            --bg-secondary: #161b22;
            --text-primary: #c9d1d9;
            --text-secondary: #8b949e;
            --accent: #58a6ff;
            --success: #3fb950;
            --warning: #d29922;
        }}
        
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem;
        }}
        
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        
        h1 {{
            font-size: 2rem;
            margin-bottom: 0.5rem;
            background: linear-gradient(135deg, var(--accent), var(--success));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }}
        
        .meta {{
            color: var(--text-secondary);
            margin-bottom: 2rem;
        }}
        
        .card {{
            background: var(--bg-secondary);
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            border: 1px solid rgba(255,255,255,0.1);
        }}
        
        .card h2 {{
            font-size: 1.25rem;
            margin-bottom: 1rem;
            color: var(--accent);
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
        }}
        
        th, td {{
            text-align: left;
            padding: 0.75rem 1rem;
            border-bottom: 1px solid rgba(255,255,255,0.1);
        }}
        
        th {{
            font-weight: 600;
            color: var(--text-secondary);
            text-transform: uppercase;
            font-size: 0.75rem;
            letter-spacing: 0.05em;
        }}
        
        .score {{
            font-family: 'SF Mono', 'Monaco', monospace;
            font-weight: 600;
        }}
        
        .score.high {{
            color: var(--success);
        }}
        
        .score.medium {{
            color: var(--warning);
        }}
        
        .score.low {{
            color: var(--text-secondary);
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
        }}
        
        .stat {{
            padding: 1rem;
            background: rgba(255,255,255,0.05);
            border-radius: 8px;
        }}
        
        .stat-value {{
            font-size: 1.5rem;
            font-weight: 700;
            color: var(--accent);
        }}
        
        .stat-label {{
            font-size: 0.875rem;
            color: var(--text-secondary);
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 PSP Benchmark Leaderboard</h1>
        <p class="meta">Generated: {leaderboard["generated_at"]}</p>
"""

    # Method cards
    for method_name, stats in sorted(leaderboard.get("methods", {}).items()):
        html += f"""
        <div class="card">
            <h2>{method_name}</h2>
            <div class="stats-grid">
                <div class="stat">
                    <div class="stat-value">{stats["total_targets"]}</div>
                    <div class="stat-label">Targets Evaluated</div>
                </div>
"""

        for scorer, scorer_stats in sorted(stats.get("scorers", {}).items()):
            mean = scorer_stats.get("mean", 0)
            score_class = "high" if mean > 0.7 else "medium" if mean > 0.5 else "low"

            html += f"""
                <div class="stat">
                    <div class="stat-value score {score_class}">{mean:.3f}</div>
                    <div class="stat-label">{scorer.upper()} (mean)</div>
                </div>
"""

        html += """
            </div>
            <table>
                <thead>
                    <tr>
                        <th>Metric</th>
                        <th>Mean</th>
                        <th>Median</th>
                        <th>Min</th>
                        <th>Max</th>
                        <th>Success Rate (>0.5)</th>
                    </tr>
                </thead>
                <tbody>
"""

        for scorer, scorer_stats in sorted(stats.get("scorers", {}).items()):
            mean = scorer_stats.get("mean", 0)
            score_class = "high" if mean > 0.7 else "medium" if mean > 0.5 else "low"

            html += f"""
                    <tr>
                        <td>{scorer.upper()}</td>
                        <td class="score {score_class}">{mean:.4f}</td>
                        <td class="score">{scorer_stats.get("median", 0):.4f}</td>
                        <td class="score">{scorer_stats.get("min", 0):.4f}</td>
                        <td class="score">{scorer_stats.get("max", 0):.4f}</td>
                        <td class="score">{scorer_stats.get("success_rate", 0) * 100:.1f}%</td>
                    </tr>
"""

        html += """
                </tbody>
            </table>
        </div>
"""

    html += """
    </div>
</body>
</html>
"""

    with open(output_path, "w") as f:
        f.write(html)

    logger.info(f"Generated HTML: {output_path}")


def _load_config(config_path):
    """Load config with local overrides."""
    import yaml
    from pathlib import Path
    config_path = Path(config_path)
    with open(config_path) as f:
        config = yaml.safe_load(f)
    local_path = config_path.parent / "config.local.yaml"
    if local_path.exists():
        with open(local_path) as f:
            local_config = yaml.safe_load(f) or {}
        for key, value in local_config.items():
            if key in config and isinstance(config[key], dict) and isinstance(value, dict):
                config[key].update(value)
            else:
                config[key] = value
    return config


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate evaluation results into leaderboard"
    )
    parser.add_argument(
        "--eval-dir",
        type=Path,
        required=True,
        help="Directory containing evaluation JSON files",
    )
    parser.add_argument(
        "--output-csv", type=Path, required=True, help="Output CSV file path"
    )
    parser.add_argument(
        "--output-html", type=Path, required=True, help="Output HTML file path"
    )
    parser.add_argument(
        "--output-json", type=Path, required=True, help="Output JSON summary file path"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="Path to configuration file",
    )

    args = parser.parse_args()

    # Load config
    config = _load_config(args.config)

    # Load results
    results = load_eval_results(args.eval_dir)

    if not results:
        logger.warning("No evaluation results found")

    # Build leaderboard
    leaderboard = build_leaderboard(results, config)

    # Ensure output directories exist
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)

    # Generate outputs
    generate_csv(results, args.output_csv)
    generate_html(leaderboard, args.output_html)

    with open(args.output_json, "w") as f:
        json.dump(leaderboard, f, indent=2)
    logger.info(f"Generated JSON: {args.output_json}")

    # Print summary
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)

    for method, stats in sorted(leaderboard.get("methods", {}).items()):
        print(f"\n{method}:")
        for scorer, scorer_stats in sorted(stats.get("scorers", {}).items()):
            print(
                f"  {scorer:12s}: mean={scorer_stats['mean']:.4f}, "
                f"success_rate={scorer_stats['success_rate'] * 100:.1f}%"
            )

    print("\n" + "=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
