import json
from src.config import Config
from src.utils import display_full_info, get_full_info, parse_args
from contextlib import redirect_stdout

if __name__ == "__main__":

    cfg = Config()

    args = parse_args()
    
    all_json_data = []

    for pdb_id in args.pdb_ids:
        full_info = get_full_info(pdb_id, cfg)
        all_json_data.append(full_info)

        if args.save_text:
            with open(args.save_text, 'a') as f:
                with redirect_stdout(f):
                    display_full_info(full_info)
        
        if args.save_json:
            with open(args.save_json, 'w') as f:
                json.dump(all_json_data, f, indent=2)

        if not args.quiet:
            display_full_info(full_info)
                        

